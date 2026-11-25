using BenchmarkTools
using OpticalSpectrumModels
using ComponentArrays
using ForwardDiff
using ExoMol

# Evaluation counts, written by tune.jl.
const PARAMS_PATH = joinpath(@__DIR__, "params.json")

const SUITE = BenchmarkGroup()

# EmissionSpectrum with 200 synthetic transitions on 500 points, also under
# ForwardDiff.
SUITE["emission_spectrum"] = let g = BenchmarkGroup()
    xs = collect(range(12400.0, 12600.0, length = 500))
    transitions = [ElectricDipoleTransition(State(:X, 0, i), State(:A, 1, i), 12500.0 + i, 1e6)
                   for i in 1:200]
    params = ComponentVector(T = 500.0)
    sp = EmissionSpectrum(;
        transitions,
        distribution_function = VariableDistribution(:T),
        lineshape_function = (t, p) -> Lorentzian(0.1),
    )
    ret = zeros(length(xs))

    g["allocating"] = @benchmarkable $sp($xs, $params)
    g["inplace"] = @benchmarkable $sp($ret, $xs, $params)

    dual_params = ComponentVector(T = ForwardDiff.Dual(500.0, 1.0))
    dual_ret = zeros(eltype(dual_params), length(xs))
    g["forwarddiff_alloc"] = @benchmarkable $sp($xs, $dual_params)
    g["forwarddiff_inplace"] = @benchmarkable $sp($dual_ret, $xs, $dual_params)
    g
end

# The case of examples/sio_emission.jl: about 255k SiO lines on 24001 points
# with LDM, with and without the slit.
SUITE["sio_emission"] = let g = BenchmarkGroup()
    SiO = load_isotopologue("SiO", "28Si-16O", "EBJT"; wn_range = (1080.0, 1320.0))
    tl = make_transition_list(SiO)
    bd = Boltzmann(SiO, upper_manifold(tl))

    m_SiO = molecular_mass_kg(SiO)
    kB = 1.380649e-23
    function sio_lineshape(t, p)
        γ = gamma_pressure(p.P / 101325.0, p.T, p.x; gamma_air = 0.07, n_air = 0.5)
        return Voigt(sigma_doppler(t.wavenumber, p.T; m = m_SiO), γ)
    end

    sp = RadianceSpectrum(;
        species = SpeciesModel(;
            transitions = tl,
            distribution_function = bd,
            lineshape_function = sio_lineshape,
            number_density_function = p -> p.P * p.x / kB / p.T / 1e6,
            method = LDM()),
        path_length_function = p -> p.L,
    )

    xs = range(1080.0, 1320.0; step = 0.01)
    params = ComponentVector(T = 1000.0, P = 1.01325e5, x = 0.1, L = 1.0)

    g["radiance"] = @benchmarkable $sp($xs, $params) evals=1 samples=5 seconds=30
    g["radiance_slit"] = @benchmarkable convolve_spectrum($sp, $xs, $params, Triangular(1.0)) evals=1 samples=5 seconds=30
    g
end

if isfile(PARAMS_PATH)
    loadparams!(SUITE, first(BenchmarkTools.load(PARAMS_PATH)), :evals)
end
