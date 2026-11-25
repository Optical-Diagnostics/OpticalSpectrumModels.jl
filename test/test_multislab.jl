using OpticalSpectrumModels
using Test
using ComponentArrays

include("synthetic.jl")

Q(T) = 1.0 + 3.0 * exp(-c2 * 1000.0 / T)
lsf(t, p) = Lorentzian(0.5)
bdk(Tkey) = Boltzmann(Q, 0.0, State(:X), p -> p[Val(Tkey)])

# One species per slab; per-slab temperature and density via captured keys.
mk_species(Tkey, Nkey; df = bdk(Tkey), lineshape_function = lsf, method = PointSample()) =
    SpeciesModel(; transitions = [t1], distribution_function = df, lineshape_function,
                   number_density_function = p -> p[Val(Nkey)], method)
slab(s, Lkey = :L) = Slab(; species = s, path_length_function = p -> p[Lkey])

@testset "Multi-slab line of sight" begin
    xs = collect(range(990.0, 1010.0; length = 401))

    @testset "Boltzmann temperature_function" begin
        p = ComponentVector(T = 800.0, TA = 800.0, TB = 1500.0)
        @test bdk(:TA)(upper0, p) == Boltzmann(Q, 0.0, State(:X))(upper0, p)   # default key :T

        # Fast-path guard: equal temperature values -> common T; unequal -> nothing.
        sA = SpeciesModel(; transitions = [t1], distribution_function = bdk(:TA), lineshape_function = lsf)
        sB = SpeciesModel(; transitions = [t1], distribution_function = bdk(:TB), lineshape_function = lsf)
        @test OpticalSpectrumModels._lte_temperature((sA, sA), p) == 800.0
        @test OpticalSpectrumModels._lte_temperature((sA, sB), p) === nothing
    end

    @testset "one Slab ≡ species form; two half-slabs ≡ one full slab" begin
        params = ComponentVector(TA = 800.0, NA = 1e17, L = 10.0)
        s = mk_species(:TA, :NA)
        r1 = RadianceSpectrum(; species = s, path_length_function = p -> p[:L])(xs, params)
        @test RadianceSpectrum(; slabs = slab(s))(xs, params) == r1
        half = Slab(; species = s, path_length_function = p -> p[:L] / 2)
        @test RadianceSpectrum(; slabs = (half, half))(xs, params) ≈ r1 rtol = 1e-12
    end

    @testset "self-reversal: hot core behind cool edge, and order matters" begin
        # Hot, optically thick back slab saturates near B(T_hot); the cooler
        # front slab absorbs the line centre and re-emits at B(T_cool) < B(T_hot)
        # -> central dip below the wings. Reversed order shows no dip.
        # τ at line centre ≈ 10 (hot) and ≈ 3 (cool); both are thin at 1002 cm⁻¹.
        params = ComponentVector(T_hot = 3000.0, T_cool = 1000.0,
                                 N_hot = 3e17,   N_cool = 2.7e16, L = 10.0)
        hot  = slab(mk_species(:T_hot,  :N_hot))
        cool = slab(mk_species(:T_cool, :N_cool))
        r_rev = RadianceSpectrum(; slabs = (hot, cool))(xs, params)   # cool in front
        r_ok  = RadianceSpectrum(; slabs = (cool, hot))(xs, params)   # hot in front
        i_c = argmin(abs.(xs .- 1000.0))          # line centre
        i_w = argmin(abs.(xs .- 1002.0))          # wing, still on the line
        @test r_rev[i_w] > r_rev[i_c]
        @test r_ok[i_c] > r_ok[i_w]
    end

    @testset "per-slab temperatures: fast path ≡ general path" begin
        # Narrow lines so the grid-ν vs rest-ν̃ Planck gap is negligible.
        mk(Tkey, Nkey, df) = slab(mk_species(Tkey, Nkey; df, lineshape_function = (t, p) -> Lorentzian(0.01),
                                             method = PointSample(wing_cut = 100.0)))
        params = ComponentVector(TA = 800.0, TB = 2000.0, NA = 1e16, NB = 5e15, L = 10.0)
        xsn = collect(range(999.5, 1000.5; length = 101))
        fast = RadianceSpectrum(; slabs = (mk(:TA, :NA, bdk(:TA)), mk(:TB, :NB, bdk(:TB))))(xsn, params)
        gen  = RadianceSpectrum(; slabs = (mk(:TA, :NA, ProductDistribution(bdk(:TA))),
                                           mk(:TB, :NB, ProductDistribution(bdk(:TB)))))(xsn, params)
        @test fast ≈ gen rtol = 1e-3
        @test any(>(0), fast)
    end

    params2 = ComponentVector(TA = 800.0, TB = 1500.0, NA = 1e16, NB = 4e15, LA = 10.0, LB = 3.0)
    slabA = slab(mk_species(:TA, :NA), :LA)
    slabB = slab(mk_species(:TB, :NB), :LB)

    @testset "linear models over slabs" begin
        A_A  = AbsorbanceSpectrum(; slabs = (slabA,))(xs, params2)
        A_B  = AbsorbanceSpectrum(; slabs = (slabB,))(xs, params2)
        A_AB = AbsorbanceSpectrum(; slabs = (slabA, slabB))(xs, params2)
        @test A_AB ≈ A_A .+ A_B rtol = 1e-12
        @test TransmittanceSpectrum(; slabs = (slabA, slabB))(xs, params2) ≈ exp.(-A_AB) rtol = 1e-12
        e_A = EmissionSpectrum(; slabs = (slabA,))(xs, params2)
        e_B = EmissionSpectrum(; slabs = (slabB,))(xs, params2)
        @test EmissionSpectrum(; slabs = (slabA, slabB))(xs, params2) ≈ e_A .+ e_B rtol = 1e-12
    end

    @testset "in-place parity and buffer reuse over slabs" begin
        for model in (RadianceSpectrum(; slabs = (slabA, slabB)),
                      TransmittanceSpectrum(; slabs = (slabA, slabB)))
            ref = model(xs, params2)
            ret = similar(ref)
            model(ret, xs, params2)
            model(ret, xs, params2)
            @test ret == ref
            a_inplace = @allocated model(ret, xs, params2)
            @test a_inplace < @allocated model(xs, params2)
        end
    end

    @testset "constructor validation" begin
        s = mk_species(:TA, :NA)
        sl = slab(s)
        @test_throws ArgumentError RadianceSpectrum(; slabs = (sl,), species = s)
        @test_throws ArgumentError RadianceSpectrum(; slabs = (sl,), path_length_function = p -> p[:L])
        @test_throws ArgumentError RadianceSpectrum(; slabs = ())
        @test_throws ArgumentError RadianceSpectrum(; species = s)   # missing path length
    end
end
