using OpticalSpectrumModels
using ExoMol
using Test
using ComponentArrays
using ForwardDiff

# Downloads to Julia scratch space on first run; cached on subsequent runs.
# Neither N2 transition file names its range, so wn_range keeps all ~8.4M lines.
const _N2  = load_isotopologue("N2", "14N2", "WCCRMT"; wn_range = (25000.0, 40000.0))
const _SiO = load_isotopologue("SiO", "28Si-16O", "EBJT"; wn_range = (1080.0, 1320.0))
const ext  = Base.get_extension(OpticalSpectrumModels, :OpticalSpectrumModelsExoMolExt)
const c2   = 100 * OpticalSpectrumModels.h * OpticalSpectrumModels.c0 / OpticalSpectrumModels.kB

loaded_states(data) = State.(ext.ExoMolState.(data.states))

@testset "N2 WCCRMT" begin
    tl = make_transition_list(_N2)
    @test all(t -> t.A > 0 && t.upper.degeneracy > 0 && t.lower.degeneracy > 0, tl)
    # ν = E_upper − E_lower checks the state lookup. The dataset's ERJ file
    # gives 810 lines a wavenumber off the states file, so check the share.
    off = count(t -> !isapprox(t.wavenumber, t.upper.energy - t.lower.energy; atol = 1e-3), tl)
    @test off / length(tl) < 1e-3

    upper_el = first(tl).upper.electronic_state
    bd = Boltzmann(_N2, State(upper_el); Ts = 200:100:5000)

    @testset "Boltzmann" begin
        # The parent is restricted, so Q is the sum over its states and the
        # populations sum to 1.
        states = loaded_states(_N2)
        @test all(T -> sum(s -> bd(s, ComponentVector(T = T)), states) ≈ 1, (300.0, 3000.0))
        @test first(tl).lower.electronic_state != upper_el
        @test bd(first(tl).lower, ComponentVector(T = 300.0)) == 0.0
    end

    @testset "EmissionSpectrum against a hand sum" begin
        sub = filter(t -> 30000 < t.wavenumber < 30050 && t.upper.electronic_state == upper_el, tl)
        @test !isempty(sub)
        xs     = collect(range(30000.0, 30050.0; length = 2001))
        params = ComponentVector(T = 3000.0)
        shape(t, T) = Voigt(sigma_doppler(t.wavenumber, T; m = 28.0 * 1.66054e-27), 0.05)
        sp = EmissionSpectrum(; transitions = sub, distribution_function = bd,
                                lineshape_function = (t, p) -> shape(t, p[:T]))
        manual = zeros(length(xs))
        for t in sub
            manual .+= t.A * bd(t.upper, params) .* shape(t, 3000.0).(xs .- t.wavenumber)
        end
        @test sp(xs, params) ≈ manual rtol = 1e-12
    end

    @testset "upper_manifold" begin
        one_manifold = filter(t -> electronic_state(t.upper) == upper_el, tl)
        @test upper_manifold(one_manifold) == State(upper_el)
        @test_throws ArgumentError upper_manifold(tl)
        @test_throws ArgumentError upper_manifold(eltype(tl)[])
    end
end

@testset "partition-function sourcing" begin
    # SiO/EBJT is the unrestricted case: a single electronic state, so
    # `State(:X)` selects every loaded state and the dataset's own total Q(T)
    # is the correct denominator.
    sio_st = loaded_states(_SiO)
    parent = State(first(sio_st).electronic_state)
    @test all(s -> s in parent, sio_st)
    bd = Boltzmann(_SiO, parent)
    # The .pf branch: Q is the file's interpolation, and E0 is the .pf file's
    # own energy zero, not the lowest included state.
    @test bd.partition_function === _SiO.partition_function
    @test bd.E0 == 0
    # The state sum misses states, so it is a lower bound on the tabulated
    # total.
    @test all((300.0, 1000.0, 3000.0)) do T
        Q_pf  = bd.partition_function(T)
        Q_sum = sum(s -> s.degeneracy * exp(-s.energy * c2 / T), sio_st)
        Q_sum <= Q_pf * (1 + 1e-9) && isapprox(Q_sum, Q_pf; rtol = 0.05)
    end

    # N2/WCCRMT holds several electronic states, so a single-state parent is
    # restricted: Q is a subset-sum spline, not the dataset's .pf.
    parent_N2 = State(first(loaded_states(_N2)).electronic_state)
    @test _N2.partition_function !== nothing
    @test Boltzmann(_N2, parent_N2; Ts = 200:100:5000).partition_function !== _N2.partition_function

    # Restricting SiO to one vibrational level sums over its states and
    # shrinks Q.
    v0 = first(s.vibrational_state for s in sio_st if s.vibrational_state !== nothing)
    bd_v0 = Boltzmann(_SiO, State(parent.electronic_state, v0))
    @test bd_v0.partition_function(1000.0) < _SiO.partition_function(1000.0)
end

@testset "state-sum Q(T) at plasma temperatures" begin
    states = [State(:X, v, J, 1240.0 * v + 0.72 * J * (J + 1), 2J + 1) for v in 0:40 for J in 0:300]
    Q, E0 = ext.partition_function(states, State(:X))
    exact(T) = sum(s.degeneracy * exp(-(s.energy - E0) * c2 / T) for s in states)
    @test all(T -> isapprox(Q(T), exact(T); rtol = 1e-6), (300.0, 1500.0, 3000.0, 5000.0))
    # Outside the grid the sum is evaluated directly, also under ForwardDiff:
    # a solver's trial step may leave the grid.
    @test all(T -> isapprox(Q(T), exact(T); rtol = 1e-12), (10.0, 20_000.0))
    e = [(s.energy - E0) * c2 for s in states]
    dQ_exact = sum(s.degeneracy * ei / 20_000.0^2 * exp(-ei / 20_000.0) for (s, ei) in zip(states, e))
    @test ForwardDiff.derivative(Q, 20_000.0) ≈ dQ_exact rtol = 1e-10
end

@testset "molecular mass from the dataset" begin
    # Read from the .def file, not summed from mass numbers.
    @test molecular_mass_kg(_SiO) ≈ 43.971842 * 1.66053906660e-27 rtol = 1e-9
    @test sigma_doppler(1252.68, 1000.0, _SiO) ==
          sigma_doppler(1252.68, 1000.0; m = molecular_mass_kg(_SiO))
end
