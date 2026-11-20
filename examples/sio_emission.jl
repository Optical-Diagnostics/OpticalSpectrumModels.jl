using SpectraUtils
using SpectraUtils: kB
using ExoMol
using ComponentArrays
using CairoMakie
using DelimitedFiles

# Reproduce the following RADIS calculation:
#
#   from radis import calc_spectrum
#   s = calc_spectrum(
#       1080, 1320,          # cm⁻¹
#       species="SiO", isotope="1",
#       pressure=1.01325,    # bar
#       Tgas=1000,           # K
#       mole_fraction=0.1, path_length=1,  # cm
#       databank=("exomol", "EBJT"),
#   )
#   s.apply_slit(1, "cm-1")
#   s.plot("radiance")
#
# RADIS "isotope 1" = ²⁸Si¹⁶O.  The window 1080–1320 cm⁻¹ covers the
# fundamental vibrational band (v=1→0, ω_e ≈ 1242 cm⁻¹).

const P      = 1.01325e5   # Pa (1.01325 bar)
const x_mol  = 0.1
const T_gas  = 1000.0      # K
const L_cm   = 1.0         # cm

const _SiO = load_isotopologue("SiO", "28Si-16O", "EBJT"; wn_range=(1080.0, 1320.0))

tl = make_transition_list(_SiO)
@info "SiO EBJT transitions in range" length(tl)
@info "Electronic states (upper)" collect(unique(electronic_state(t.upper) for t in tl))

# At 1000 K, SiO occupies only the ground electronic state X¹Σ.
const parent = upper_manifold(tl)
const bd     = Boltzmann(_SiO, parent)   # SiO/EBJT is single-electronic-state, so Q(T) comes from its .pf file

# ²⁸Si¹⁶O molecular mass from the dataset, in kg
const m_SiO = molecular_mass_kg(_SiO)

# Pressure broadening. SiO EBJT has no per-line `.broad` file, so RADIS and
# this example use the dataset-wide default from its `.def.json`, for both air
# and self broadening.
const γ_ref  = 0.07   # cm⁻¹/atm at Tref, SiO EBJT def.json default (air = self)
const n_exp  = 0.5    # temperature exponent, same source
const Tref_K = 296.0  # K, ExoMol/HITRAN reference temperature
const atm_Pa = 101325.0  # Pa per atm

function sio_lineshape(t, p)
    T_, P_, x_ = p.T, p.P, p.x
    p_atm = P_ / atm_Pa
    γ = gamma_pressure(p_atm * (1 - x_), T_; gamma_ref=γ_ref, n=n_exp, Tref=Tref_K) +
        gamma_pressure(p_atm * x_,       T_; gamma_ref=γ_ref, n=n_exp, Tref=Tref_K)
    return Voigt(sigma_doppler(t.wavenumber, T_; m=m_SiO), γ)
end

# RADIS's `apply_slit(1, "cm-1")` is triangular. A triangular slit has no
# closed form with Voigt, so `instrument_function` stays unset and the slit
# is applied after the radiative transfer with `convolve_spectrum`.
#
# `LDM()` matches RADIS's default here. `PointSample()` without a `wing_cut`
# evaluates all 254,675 lines at all 24,001 grid points, about 50× slower.
sp = RadianceSpectrum(;
    species = SpeciesModel(;
        transitions              = tl,
        distribution_function    = bd,
        lineshape_function       = sio_lineshape,
        number_density_function  = (p) -> p.P * p.x / kB / p.T / 1e6,  # Pa·mol⁻¹·K⁻¹ → cm⁻³
        method                   = LDM()),
    path_length_function = (p) -> p.L,
)

# 0.01 cm⁻¹ grid, RADIS's default `wstep` here. A finer grid (e.g. 0.005 cm⁻¹)
# resolves the Voigt linewidth (~0.077 cm⁻¹ FWHM) better.
xs     = range(1080.0, 1320.0; step=0.01)
params = ComponentVector(T=T_gas, P=P, x=x_mol, L=L_cm)

radiance_W = convolve_spectrum(sp, xs, params, Triangular(1.0))  # W cm⁻² sr⁻¹ (cm⁻¹)⁻¹
radiance   = radiance_W .* 1e3  # mW cm⁻² sr⁻¹ (cm⁻¹)⁻¹
@info "Radiance computed" length=length(xs) peak_mW_per_cm2_sr_wn=maximum(radiance)

csv_path = joinpath(@__DIR__, "sio_emission.csv")
open(csv_path, "w") do io
    writedlm(io, ["wavenumber_cm1" "radiance_mW_cm2_sr_cm1"], ',')
    writedlm(io, [xs radiance], ',')
end
@info "CSV written" path=csv_path rows=length(xs)

pdf_path = joinpath(@__DIR__, "sio_emission.pdf")
fig = CairoMakie.Figure(size=(800, 400))
ax  = CairoMakie.Axis(fig[1, 1];
    xlabel = "Wavenumber [cm⁻¹]",
    ylabel = "Radiance [mW cm⁻² sr⁻¹ (cm⁻¹)⁻¹]",
    title  = "SiO emission, ExoMol EBJT, T = $(Int(T_gas)) K, P = 1.01325 bar, x = $x_mol",
)
CairoMakie.lines!(ax, xs, radiance; linewidth=0.8)
CairoMakie.save(pdf_path, fig)
@info "PDF written" path=pdf_path
