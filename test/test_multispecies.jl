using OpticalSpectrumModels
using Test
using ComponentArrays

include("synthetic.jl")

# Two species, one carrying t1 and one t2 (overlapping wings), with
# independent number densities in params.
Q1(T) = 1.0 + 3.0 * exp(-c2 * 1000.0 / T)
Q2(T) = 1.0 + 5.0 * exp(-c2 * 1050.0 / T)
const bd1 = Boltzmann(Q1, 0.0, State(:X))
const bd2 = Boltzmann(Q2, 0.0, State(:X))

lsf(t, p) = Lorentzian(0.5)

mk_species(tr, df, key; method = PointSample(), lineshape_function = lsf) =
    SpeciesModel(; transitions = [tr], distribution_function = df, lineshape_function,
                   number_density_function = p -> p[Val(key)], method)

@testset "Multi-species in one slab" begin
    xs     = collect(range(950.0, 1100.0; length = 300))
    params = ComponentVector(T = 800.0, N1 = 2.0e15, N2 = 5.0e14, L = 10.0)
    PL     = (path_length_function = p -> p[:L],)

    s1 = mk_species(t1, bd1, :N1)
    s2 = mk_species(t2, bd2, :N2)
    α1 = AbsorptionCoefficientSpectrum(; species = s1)(xs, params)
    α2 = AbsorptionCoefficientSpectrum(; species = s2)(xs, params)

    @testset "linear models add" begin
        @test AbsorptionCoefficientSpectrum(; species = (s1, s2))(xs, params) ≈ α1 .+ α2 rtol = 1e-12
        @test TransmittanceSpectrum(; species = (s1, s2), PL...)(xs, params) ≈
              exp.(-(α1 .+ α2) .* 10.0) rtol = 1e-12
        e1 = EmissionSpectrum(; species = s1, PL...)(xs, params)
        e2 = EmissionSpectrum(; species = s2, PL...)(xs, params)
        @test EmissionSpectrum(; species = (s1, s2), PL...)(xs, params) ≈ e1 .+ e2 rtol = 1e-12
    end

    @testset "radiance: LTE fast path over combined α" begin
        r12 = RadianceSpectrum(; species = (s1, s2), PL...)(xs, params)
        @test r12 ≈ @.(planck_radiance(xs, 800.0) * -expm1(-(α1 + α2) * 10.0)) rtol = 1e-12
    end

    @testset "radiance: mixed LTE/non-LTE takes the general path and agrees" begin
        # Wrapping one species' Boltzmann forces the general (α, j) path for the
        # whole slab; narrow lines make the two paths agree (see test_radiance.jl).
        nlsf(t, p) = Lorentzian(0.01)
        n1 = mk_species(t1, bd1, :N1; lineshape_function = nlsf)
        n2(df) = mk_species(t2, df, :N2; lineshape_function = nlsf)
        xsn = sort(vcat(range(999.8, 1000.2; length = 41), range(1049.8, 1050.2; length = 41)))
        r_lte   = RadianceSpectrum(; species = (n1, n2(bd2)), PL...)(xsn, params)
        r_mixed = RadianceSpectrum(; species = (n1, n2(ProductDistribution(bd2))), PL...)(xsn, params)
        @test r_mixed ≈ r_lte rtol = 1e-3
    end

    @testset "combined radiance is not the sum of single-species radiances" begin
        # Two species sharing one optically thick line: the slab saturates at
        # B(T), while a sum of radiances would double it.
        thick = ComponentVector(T = 800.0, N1 = 1e19, N2 = 0.0, L = 100.0)
        r_comb = RadianceSpectrum(; species = (s1, mk_species(t1, bd1, :N1)), PL...)([1000.0], thick)
        r_one  = RadianceSpectrum(; species = s1, PL...)([1000.0], thick)
        @test r_comb[1] ≈ r_one[1] rtol = 1e-10
    end

    @testset "per-species methods coexist" begin
        # Δx = 0.1 × hwhm, so BinAverage's bin average matches PointSample's point value.
        xsf = collect(range(950.0, 1100.0; length = 3001))
        s1m = mk_species(t1, bd1, :N1; method = PointSample(wing_cut = 100.0))
        s2m = mk_species(t2, bd2, :N2; method = BinAverage(bin_width = xsf[2] - xsf[1]))
        α_mix = AbsorptionCoefficientSpectrum(; species = (s1m, s2m))(xsf, params)
        α_ref = AbsorptionCoefficientSpectrum(; species = (s1, s2))(xsf, params)
        @test maximum(abs, α_mix .- α_ref) / maximum(α_ref) < 1e-3
    end

    @testset "in-place parity" begin
        for model in (AbsorptionCoefficientSpectrum(; species = (s1, s2)),
                      RadianceSpectrum(; species = (s1, s2), PL...))
            ref = model(xs, params)
            ret = similar(ref)
            model(ret, xs, params)
            model(ret, xs, params)
            @test ret == ref
        end
    end

    @testset "constructor validation" begin
        @test_throws ArgumentError EmissionSpectrum(; species = s1, transitions = [t1])
        @test_throws ArgumentError EmissionSpectrum(; species = ())
        # `method` belongs to SpeciesModel; the flat keywords use PointSample().
        @test_throws MethodError EmissionSpectrum(;
            transitions = [t1], distribution_function = bd1,
            lineshape_function = lsf, method = PointSample(wing_cut = 50.0))
        @test EmissionSpectrum(; transitions = [t1], distribution_function = bd1,
            lineshape_function = lsf).slabs[1].species[1].method == PointSample()
    end
end
