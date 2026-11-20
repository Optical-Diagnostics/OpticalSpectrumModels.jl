using SpectraUtils
using CairoMakie

# Transmittance of a gas layer from a HITRAN or HITEMP `.par` file.
#
#     julia --project=examples examples/hitran_transmittance.jl [file.par]
#
# Download the lines from https://hitran.org (".par" output format), for
# example CO (molecule 5) between 2000 and 2300 cm⁻¹. Without an argument the
# example reads the four-line CO test fixture.
#
# HITRAN's line intensity S already includes the natural isotopic abundance,
# so each line's absorption coefficient is S(T)·N·φ(ν), with N the number
# density of the molecule.
const par_file = isempty(ARGS) ? joinpath(@__DIR__, "..", "test", "fixtures", "CO_sample.par") : ARGS[1]

const T     = 1000.0      # K
const P_atm = 1.0         # atm
const x     = 0.05        # mole fraction
const L     = 10.0        # cm
const N     = P_atm * 101325 * x / (SpectraUtils.kB * T) / 1e6   # molecules/cm³
const Tref  = 296.0       # K, HITRAN reference temperature

par_lines = HITRAN.read_par(par_file; wn_min = 2000.0, wn_max = 2300.0)
@info "Lines read" file = par_file n = length(par_lines)

# Air and self broadening from the line list. HITRAN gives one temperature
# exponent, for air; it is used for both here.
function absorption_line(l)
    S = HITRAN.linestrength(l, T; Tref)
    γ = gamma_pressure(P_atm, T, x; gamma_air = l.airbrd, n_air = l.Tdpair,
                                gamma_self = l.selbrd, Tref)
    σ = sigma_doppler(l.wav, T, HITRAN.molparam(l.mol_id, l.iso))
    return Line(l.wav, S * N, Voigt(σ, γ))
end
ls = map(absorption_line, par_lines)

xs = range(2000.0, 2300.0; step = 0.005)
α  = calculate_lines(xs, ls, SparseLDM())          # cm⁻¹
transmittance = exp.(-α .* L)
@info "Transmittance range" min = minimum(transmittance) max = maximum(transmittance)

pdf_path = joinpath(@__DIR__, "hitran_transmittance.pdf")
fig = Figure(size = (1000, 400))
ax  = Axis(fig[1, 1]; xlabel = "Wavenumber [cm⁻¹]", ylabel = "Transmittance",
           title = "$(basename(par_file)): $(Int(T)) K, $(P_atm) atm, x = $x, L = $L cm")
lines!(ax, xs, transmittance; linewidth = 0.5)
save(pdf_path, fig)
@info "PDF written" path = pdf_path
