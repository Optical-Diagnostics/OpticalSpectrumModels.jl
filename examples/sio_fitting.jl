using OpticalSpectrumModels
using ExoMol
using ComponentArrays
using FittingUtils
using NonlinearSolve
using Distributions
using Random
using CairoMakie

# Fits gas temperature and SiO mole fraction to a synthetic radiance spectrum,
# with the same ExoMol SiO EBJT physics as `sio_emission.jl`. Restricted to a
# 1085-1095 cm⁻¹ sub-window (930 of the band's 254,675 transitions):
# `LevenbergMarquardtFit` needs a `ForwardDiff`-differentiable model. All four
# synthesis methods are; `PointSample()` is the cheapest at this line count.

Random.seed!(42)

const P_fixed = 1.01325e5  # Pa (1.01325 bar), held fixed
const L_fixed = 1.0        # cm, held fixed

const _SiO = load_isotopologue("SiO", "28Si-16O", "EBJT"; wn_range=(1080.0, 1320.0))
tl_all = make_transition_list(_SiO)
tl = filter(t -> 1085.0 <= t.wavenumber <= 1095.0, tl_all)
@info "SiO transitions in fitting window" length(tl)

const parent = upper_manifold(tl)
const bd     = Boltzmann(_SiO, parent)   # SiO/EBJT is single-electronic-state, so Q(T) comes from its .pf file
const m_SiO  = molecular_mass_kg(_SiO)  # kg

# Same broadening convention as `sio_emission.jl`: SiO EBJT has no per-line
# `.broad` file, so RADIS and this example use its dataset-wide default.
const γ_ref, n_exp, Tref_K, atm_Pa = 0.07, 0.5, 296.0, 101325.0

function sio_lineshape(t, p)
    p_atm = p.P / atm_Pa
    γ = gamma_pressure(p_atm * (1 - p.x), p.T; gamma_ref=γ_ref, n=n_exp, Tref=Tref_K) +
        gamma_pressure(p_atm * p.x,       p.T; gamma_ref=γ_ref, n=n_exp, Tref=Tref_K)
    return Voigt(sigma_doppler(t.wavenumber, p.T; m=m_SiO), γ)
end

sp = RadianceSpectrum(;
    species = SpeciesModel(;
        transitions              = tl,
        distribution_function    = bd,
        lineshape_function       = sio_lineshape,
        number_density_function  = (p) -> p.P * p.x / OpticalSpectrumModels.kB / p.T / 1e6),  # Pa·mol⁻¹·K⁻¹ → cm⁻³
    path_length_function = (p) -> p.L,
)

xs = range(1085.0, 1095.0; step=0.005)

# `fit_model` passes only the fitted parameters (T, x) to the model, so P and
# L are closed over from the constants above.
function sio_spectrum(xs, params)
    p = ComponentVector(T=params.T, P=P_fixed, x=params.x, L=L_fixed)
    return sp(xs, p) .* 1e3  # W → mW cm⁻² sr⁻¹ (cm⁻¹)⁻¹
end

# Synthetic data
const true_params = ComponentVector(T=1000.0, x=0.1)
true_spectrum = sio_spectrum(xs, true_params)

const σ_noise = 0.05 * maximum(true_spectrum)
noisy_spectrum = true_spectrum .+ rand(Normal(0, σ_noise), length(xs))
data = Normal.(noisy_spectrum, σ_noise)

@info "Synthetic data generated" n_points=length(xs) σ_noise peak=maximum(true_spectrum)

# Priors double as initial guesses (their means) for LevenbergMarquardtFit.
priors = ComponentVector(T=Normal(900.0, 200.0), x=Normal(0.15, 0.05))

res = fit_model(sio_spectrum, xs, data, priors, LevenbergMarquardtFit())

@info "Fit results" T=coef(res).T T_stderr=stderror(res).T x=coef(res).x x_stderr=stderror(res).x
@info "True values" T=true_params.T x=true_params.x

best_fit = sio_spectrum(xs, coef(res))

fig = CairoMakie.Figure(size=(800, 500))
ax = CairoMakie.Axis(fig[1, 1];
    xlabel = "Wavenumber [cm⁻¹]",
    ylabel = "Radiance [mW cm⁻² sr⁻¹ (cm⁻¹)⁻¹]",
    title  = "SiO fit: T = $(round(coef(res).T, digits=1)) K (true 1000 K), " *
             "x = $(round(coef(res).x, digits=3)) (true 0.1)",
)
CairoMakie.scatter!(ax, xs, noisy_spectrum; markersize=3, color=(:gray, 0.5), label="Synthetic data")
CairoMakie.lines!(ax, xs, true_spectrum; color=:black, linestyle=:dash, linewidth=1.5, label="True")
CairoMakie.lines!(ax, xs, best_fit; color=:red, linewidth=1.5, label="Best fit")
CairoMakie.axislegend(ax)

pdf_path = joinpath(@__DIR__, "sio_fitting.pdf")
CairoMakie.save(pdf_path, fig)
@info "PDF written" path=pdf_path
