using OpticalSpectrumModels
using OpticalSpectrumModels: h, c0, kB
using ExoMol
using ComponentArrays
using DelimitedFiles
using CairoMakie

# N₂ Second Positive System: C³Πu → B³Πg (~315-434 nm, strongest band head at
# 337.2 nm). Reproduces an OES-toolbox reference spectrum (Trot=300 K,
# Tvib=1200 K, pseudo-Voigt instrument profile: shape η=0.5, width=0.01 nm),
# see examples/data/README.md and examples/data/OES-toolbox_SPS_Trot=300_Tvib=1200.csv.

const c2 = 100 * h * c0 / kB   # cm·K, second radiation constant (state energies are in cm⁻¹)

# N2/WCCRMT ships its transitions across two files: the main
# `14N2__WCCRMT.trans.bz2` and `14N2__WCCRMT__ERJ.trans.bz2`, which holds the
# C³Πu transitions this example needs. Neither is reduced by `wn_range` (the
# main file is unsegmented, `__ERJ`'s range is not encoded in its name), so it
# only documents the intended window here.
const _N2 = load_isotopologue("N2", "14N2", "WCCRMT"; wn_range=(20000.0, 35000.0))

tl_all = make_transition_list(_N2)
@info "Total N2/WCCRMT transitions" length(tl_all)

# Second positive system: C³Πu (upper) → B³Πg (lower). State labels from
# WCCRMT: :C3Piu / :B3Pig. Restrict to the band's actual extent (315-434 nm →
# ~23000-31700 cm⁻¹) with margin.
tl = filter(tl_all) do t
    electronic_state(t.upper) == Symbol("C3Piu") &&
    electronic_state(t.lower) == Symbol("B3Pig") &&
    20000.0 <= t.wavenumber <= 35000.0
end
@info "C-B (Second Positive System) transitions in range" length(tl)

# Reposition each transition in wavelength (nm). The reference axis is nm and
# its instrument profile has a fixed width in nm, which is a uniform
# convolution kernel only on an nm grid. `sigma_doppler` is a fractional width
# (f0 * sqrt(kB*T/m)/c0), so it applies to a wavelength position unchanged.
tl_nm = [ElectricDipoleTransition(t.lower, t.upper, 1e7 / t.wavenumber, t.A) for t in tl]


# Two-temperature (Trot ≠ Tvib) population
# Trot sets the rotational distribution within each vibrational level, Tvib
# the distribution across levels. ExoMol tabulates only the total term energy of each (v, J)
# state, so the lowest energy at each v is taken as the band origin and the
# energy above it as the rotational part.
upper_states = unique(t.upper for t in tl_nm)
vibs = sort(unique(vibrational_state(s) for s in upper_states))
states_by_v = Dict(v => filter(s -> vibrational_state(s) == v, upper_states) for v in vibs)
Evib = Dict(v => minimum(energy(s) for s in states_by_v[v]) for v in vibs)
@info "Vibrational levels (C³Πu)" vibs

struct TwoTemperatureBoltzmann{EV,SV}
    Evib::EV            # Dict{vib_type, Float64}
    states_by_v::SV      # Dict{vib_type, Vector{<:State}}
end

function (d::TwoTemperatureBoltzmann)(s, params)
    v = vibrational_state(s)
    haskey(d.Evib, v) || return zero(eltype(params))
    Tvib = params[Val(:Tvib)]
    Trot = params[Val(:Trot)]

    Qvib = sum(exp(-Ev * c2 / Tvib) for Ev in values(d.Evib))
    E0 = d.Evib[v]
    Qrot = sum(degeneracy(s2) * exp(-(energy(s2) - E0) * c2 / Trot) for s2 in d.states_by_v[v])

    return exp(-E0 * c2 / Tvib) / Qvib *
           degeneracy(s) * exp(-(energy(s) - E0) * c2 / Trot) / Qrot
end

const bd = TwoTemperatureBoltzmann(Evib, states_by_v)

# OES-toolbox's pseudo-Voigt takes a mixing fraction η (its "shape", 0 for
# Gaussian, 1 for Lorentzian) and one shared FWHM. OpticalSpectrumModels's `PseudoVoigt`
# derives η from separate Gaussian and Lorentzian widths, so a local callable
# is used here. `convolve_spectrum` needs only a callable with an `fwhm`
# method.
struct SharedFWHMPseudoVoigt
    fwhm::Float64
    eta::Float64
end
function (pv::SharedFWHMPseudoVoigt)(x)
    hwhm = pv.fwhm / 2
    sigma = pv.fwhm / (2 * sqrt(2 * log(2)))
    return pv.eta * Lorentzian(hwhm)(x) + (1 - pv.eta) * Gaussian(sigma)(x)
end
OpticalSpectrumModels.fwhm(pv::SharedFWHMPseudoVoigt) = pv.fwhm

const instrument = SharedFWHMPseudoVoigt(0.01, 0.5)  # width=0.01 nm, shape η=0.5

# The Doppler FWHM (~8e-4 nm) is below the 0.001 nm grid step, so point
# sampling would miss or overweight lines.
xs = range(300.0, 450.0; step=0.001)
sp = EmissionSpectrum(;
    species = SpeciesModel(;
        transitions           = tl_nm,
        distribution_function = bd,
        # Doppler broadening at the gas temperature, assumed equal to Trot:
        # rotation equilibrates with translation much faster than vibration.
        lineshape_function    = (t, p) -> Gaussian(sigma_doppler(t.wavenumber, p.Trot, _N2)),
        method                = BinAverage(bin_width=step(xs), wing_cut=100.0)),
)

params = ComponentVector(Trot=300.0, Tvib=1200.0)

spectrum = convolve_spectrum(sp, xs, params, instrument)
spectrum ./= maximum(spectrum)  # OES-toolbox's reference is peak-normalized
@info "Spectrum computed" length=length(xs) peak_at=xs[argmax(spectrum)]

csv_path = joinpath(@__DIR__, "sps_emission.csv")
open(csv_path, "w") do io
    writedlm(io, ["wavelength_nm" "intensity_normalized"], ',')
    writedlm(io, [xs spectrum], ',')
end
@info "CSV written" path=csv_path rows=length(xs)

pdf_path = joinpath(@__DIR__, "sps_emission.pdf")
fig = CairoMakie.Figure(size=(800, 400))
ax  = CairoMakie.Axis(fig[1, 1];
    xlabel = "Wavelength [nm]",
    ylabel = "Intensity (normalized)",
    title  = "N₂ Second Positive System (C-B), Trot = 300 K, Tvib = 1200 K",
)
CairoMakie.lines!(ax, xs, spectrum; linewidth=0.8)
CairoMakie.save(pdf_path, fig)
@info "PDF written" path=pdf_path
