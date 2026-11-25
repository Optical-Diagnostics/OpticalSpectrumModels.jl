using OpticalSpectrumModels
using OpticalSpectrumModels: HITRAN
using Test
using ComponentArrays

const CO_FIXTURE = joinpath(@__DIR__, "fixtures", "CO_sample.par")
const CO2_FIXTURE = joinpath(@__DIR__, "fixtures", "CO2_sample.par")

@testset "read_par" begin
    lines = HITRAN.read_par(CO_FIXTURE)
    @test length(lines) == 4
    @test HITRAN.molecule_name(lines[1].mol_id) == "CO"
    @test all(l -> l.mol_id == 5 && l.iso == 1, lines)

    # Every fixed-width field of the first record.
    l1 = lines[1]
    @test l1.wav ≈ 2100.0
    @test l1.intensity ≈ 1.0e-23
    @test l1.A ≈ 3.5
    @test l1.airbrd ≈ 0.070
    @test l1.selbrd ≈ 0.100
    @test l1.El ≈ 10.0
    @test l1.Tdpair ≈ 0.75
    @test l1.Pshft ≈ -0.005
    @test strip(l1.locu) == "P  1"
    @test l1.lmix == ' '
    @test l1.gp ≈ 21.0
    @test l1.gpp ≈ 19.0

    sub = HITRAN.read_par(CO_FIXTURE; wn_min = 2140.0, wn_max = 2210.0)
    @test Set(round.(Int, [l.wav for l in sub])) == Set([2150, 2200])

    # CO2 isotopologue "A" means 11.
    @test [l.iso for l in HITRAN.read_par(CO2_FIXTURE)] == [1, 11]

    @test length(HITRAN.read_hitemp([CO_FIXTURE, CO_FIXTURE])) == 8

    # A hitran.org export can mix molecules, so the id is read per record.
    mixed = joinpath(mktempdir(), "mixed.par")
    write(mixed, read(CO_FIXTURE), read(CO2_FIXTURE))
    @test [l.mol_id for l in HITRAN.read_par(mixed)] == [5, 5, 5, 5, 2, 2]
end

@testset "molparam" begin
    p = HITRAN.molparam(5, 1)   # ¹²C¹⁶O, the dominant CO isotopologue
    @test p.abundance ≈ 0.986544 rtol = 1e-9
    @test p.molar_mass ≈ 27.994915 atol = 1e-6
    @test p.gj == 1
    @test HITRAN.molparam(2, 1).isoname == "626"
end

@testset "TIPS partition sums" begin
    Q = HITRAN.partition_function(2, 1)   # CO2, 626
    @test HITRAN.temperature_range(Q) == (1.0, 5000.0)
    # Reference values from HAPI 1.3 `PYTIPS2025`: a grid point, the
    # four-point interior, and the three-point ends of the table.
    @test Q(300.0) ≈ 291.0406 rtol = 1e-7
    @test Q(296.0) ≈ 286.0939488 rtol = 1e-7
    @test Q(3.3) ≈ 3.181571785 rtol = 1e-7
    @test Q(4995.5) ≈ 884295.2129 rtol = 1e-7
    @test HITRAN.partition_function(1, 3)(296.0) ≈ 1052.144553 rtol = 1e-7
    @test HITRAN.partition_function(5, 1)(296.0) ≈ 107.4205072 rtol = 1e-7
    @test_throws DomainError Q(0.5)
    @test_throws DomainError Q(5000.1)
    @test_throws KeyError HITRAN.partition_function(2, 99)
    @test !haskey(HITRAN.TIPS, (34, 1))   # atomic O: no TIPS data

    # Consistent with the Q(296 K) column of molparam.txt for the main
    # molecules. That column comes from an earlier TIPS edition; later
    # editions revised several minor species, e.g. CH3OH by a factor of 2.
    @test all(HITRAN.MOLPARAM) do (key, p)
        key[1] <= 21 && haskey(HITRAN.TIPS, key) || return true
        isapprox(HITRAN.TIPS[key](296.0), p.Q296K; rtol = 0.02)
    end
end

@testset "linestrength" begin
    Q = HITRAN.partition_function(5, 1)
    @test HITRAN.linestrength(1.0e-23, 10.0, 2100.0, Q, 296.0) == 1.0e-23
    # Against the formula written out, with CODATA c2 = 1.438777 cm K.
    c2 = 1.4387768775039337
    S0, El, wav, T = 1.0e-23, 500.0, 2100.0, 1000.0
    manual = S0 * Q(296.0) / Q(T) * exp(-c2 * El / T) / exp(-c2 * El / 296.0) *
             (1 - exp(-c2 * wav / T)) / (1 - exp(-c2 * wav / 296.0))
    @test HITRAN.linestrength(S0, El, wav, Q, T) ≈ manual rtol = 1e-12
    l = first(HITRAN.read_par(CO_FIXTURE))
    @test HITRAN.linestrength(l, 500.0) ≈ HITRAN.linestrength(l.intensity, l.El, l.wav, Q, 500.0)
end

@testset "make_transition_list" begin
    tl = make_transition_list(HITRAN.read_par(CO_FIXTURE))
    @test length(tl) == 4
    t = first(tl)
    @test t.wavenumber ≈ 2100.0
    @test t.A ≈ 3.5
    @test t.lower.energy ≈ 10.0
    @test t.upper.energy ≈ 10.0 + 2100.0
    @test t.lower.degeneracy ≈ 19.0
    @test t.upper.degeneracy ≈ 21.0
    @test t.pressure_shift ≈ -0.005

    # One SpeciesModel applies one distribution and number density to all its
    # lines; the CO2 fixture holds isotopologues 1 and 11.
    co2 = HITRAN.read_par(CO2_FIXTURE)
    @test_logs (:warn, r"span 2 isotopologues") make_transition_list(co2)
    @test_logs make_transition_list(filter(l -> l.iso == 1, co2))

    # Unassigned lines (no A, zero weights, El = -1) would make α NaN.
    good = HITRAN.read_par(CO_FIXTURE)
    edit(l; kw...) = HITRAN.Line((get(kw, f, getfield(l, f)) for f in fieldnames(HITRAN.Line))...)
    @test all(HITRAN.is_usable, good)
    for bad in (edit(good[1]; A = NaN), edit(good[1]; A = 0.0), edit(good[1]; gp = 0.0),
                edit(good[1]; gpp = 0.0), edit(good[1]; El = -1.0))
        @test !HITRAN.is_usable(bad)
        @test_throws ArgumentError make_transition_list([good; bad])
        @test length(make_transition_list(filter(HITRAN.is_usable, [good; bad]))) == length(good)
    end
end

@testset "EmissionSpectrum end-to-end" begin
    tl = make_transition_list(HITRAN.read_par(CO_FIXTURE))
    bd = Boltzmann(HITRAN.partition_function(5, 1), 0.0, State(:X))
    xs = collect(range(2090.0, 2260.0; length = 2001))
    params = ComponentVector(T = 300.0)
    sp = EmissionSpectrum(; transitions = tl, distribution_function = bd,
                            lineshape_function = (t, p) -> Lorentzian(0.05))

    # Photon-rate amplitude A·f_u on each line.
    manual = zeros(length(xs))
    for t in tl
        manual .+= t.A * bd(t.upper, params) .* Lorentzian(0.05).(xs .- t.wavenumber)
    end
    @test sp(xs, params) ≈ manual rtol = 1e-12
end

@testset "pressure shift" begin
    tl = make_transition_list(HITRAN.read_par(CO_FIXTURE))[1:1]
    bd = Boltzmann(HITRAN.partition_function(5, 1), 0.0, State(:X))
    xs = range(2099.0, 2101.0; step = 0.001)
    shift(t, p) = t.pressure_shift * p.P
    make(lineshape) = EmissionSpectrum(; species = SpeciesModel(;
        transitions = tl, distribution_function = bd, lineshape_function = lineshape,
        lineshift_function = shift, method = PointSample()))
    peak(y) = xs[argmax(y)]

    y0 = make((t, p) -> Gaussian(0.02))(xs, ComponentVector(T = 300.0, P = 0.0))
    y1 = make((t, p) -> Gaussian(0.02))(xs, ComponentVector(T = 300.0, P = 10.0))
    @test peak(y0) ≈ 2100.0 atol = 0.0015
    @test peak(y1) ≈ 2100.0 - 0.05 atol = 0.0015
    @test sum(y1) ≈ sum(y0) rtol = 1e-6   # a shift moves the line, nothing else

    # With static lineshapes the centres must still follow the parameters.
    sp = make(StaticLineshapes(t -> Gaussian(0.02)))
    @test peak(sp(xs, ComponentVector(T = 300.0, P = 10.0))) ≈ 2100.0 - 0.05 atol = 0.0015
    @test peak(sp(xs, ComponentVector(T = 300.0, P = 0.0))) ≈ 2100.0 atol = 0.0015

    ab = AbsorbanceSpectrum(; species = SpeciesModel(; transitions = tl, distribution_function = bd,
        lineshape_function = (t, p) -> Gaussian(0.02), lineshift_function = shift, method = PointSample(),
        number_density_function = p -> 1e15), path_length_function = p -> 1.0)
    @test peak(ab(xs, ComponentVector(T = 300.0, P = 10.0))) ≈ 2100.0 - 0.05 atol = 0.0015
end

@testset "HITRAN.Line does not shadow Line" begin
    # If both exported `Line`, it would be ambiguous after
    # `using OpticalSpectrumModels.HITRAN`.
    m = Module()
    Core.eval(m, :(using OpticalSpectrumModels, OpticalSpectrumModels.HITRAN))
    @test Core.eval(m, :Line) === OpticalSpectrumModels.Line
end
