using OpticalSpectrumModels
using Test
using ComponentArrays
using ForwardDiff

include("synthetic.jl")

# The same shapes, declared both ways.
shape(t) = Voigt(2e-3 * t.wavenumber / 1000.0, 0.05)
const lsf_dyn = (t, p) -> shape(t)
const lsf_sta = StaticLineshapes(shape)

mk(lsf; kw...) = EmissionSpectrum(;
    transitions = transitions, distribution_function = bd,
    lineshape_function = lsf, kw...)

@testset "StaticLineshapes" begin
    xs = collect(range(990.0, 1060.0; length = 400))

    @testset "a params-dependent instrument is followed, also under AD" begin
        instr = (wn, p) -> Gaussian(p.w)
        sp_d = mk(lsf_dyn; instrument_function = instr)
        sp_s = mk(lsf_sta; instrument_function = instr)
        for w in (0.05, 0.5, 0.05)
            p = ComponentVector(T = 900.0, w = w)
            @test sp_s(xs, p) ≈ sp_d(xs, p) rtol = 1e-12
        end
        p = ComponentVector(T = 900.0, w = 0.2)
        @test ForwardDiff.jacobian(q -> sp_s(xs, q), p) ≈
              ForwardDiff.jacobian(q -> sp_d(xs, q), p) rtol = 1e-10
    end

    @testset "static ≡ dynamic across parameter changes" begin
        sp_d = mk(lsf_dyn)
        sp_s = mk(lsf_sta)
        for T in (800.0, 1500.0)   # repeated calls must track params
            p = ComponentVector(T = T)
            @test sp_s(xs, p) ≈ sp_d(xs, p) rtol = 1e-13
        end
        instr = (wn, p) -> Gaussian(0.02)
        p = ComponentVector(T = 900.0)
        @test mk(lsf_sta; instrument_function = instr)(xs, p) ≈
              mk(lsf_dyn; instrument_function = instr)(xs, p) rtol = 1e-13
    end

    @testset "absorption-side models" begin
        p = ComponentVector(T = 700.0, N = 1e15, L = 5.0)
        nl = (number_density_function = (p) -> p[:N],)
        a(lsf) = AbsorptionCoefficientSpectrum(; transitions, distribution_function = bd,
                                                 lineshape_function = lsf, nl...)
        @test a(lsf_sta)(xs, p) ≈ a(lsf_dyn)(xs, p) rtol = 1e-13

        # LTE fast path and the non-LTE (α, j) path.
        for df in (bd, ProductDistribution(bd))
            r(lsf) = RadianceSpectrum(; transitions, distribution_function = df,
                                        lineshape_function = lsf, nl...,
                                        path_length_function = (p) -> p[:L])
            @test r(lsf_sta)(xs, p) ≈ r(lsf_dyn)(xs, p) rtol = 1e-13
        end
    end

    @testset "gradients match; shapes stay Float64 under duals" begin
        sp_d = mk(lsf_dyn)
        sp_s = mk(lsf_sta)
        p = ComponentVector(T = 900.0)
        sp_s(xs, p)   # warm the Float64 buffer
        @test ForwardDiff.gradient(pp -> sum(sp_s(xs, pp)), p) ≈
              ForwardDiff.gradient(pp -> sum(sp_d(xs, pp)), p) rtol = 1e-12

        # The stored static lines under a Dual pass: Float64 position/shape,
        # Dual amplitude.
        pd = ComponentVector(T = ForwardDiff.Dual(900.0, 1.0))
        ls = OpticalSpectrumModels._fill_species_emission!(
            sp_s.species_caches[1][1], sp_s.slabs[1].species[1], sp_s.instrument_function, sp_s.amplitude, pd)
        @test eltype(ls).parameters[1] == Float64                     # Tx (position/shape)
        @test eltype(ls).parameters[2] <: ForwardDiff.Dual            # Ta (amplitude)
    end

    @testset "update path reuses the buffer" begin
        sp_s = mk(lsf_sta)
        s = sp_s.slabs[1].species[1]
        c = sp_s.species_caches[1][1]
        refill(T) = OpticalSpectrumModels._fill_species_emission!(c, s, sp_s.instrument_function, sp_s.amplitude,
                                                       ComponentVector(T = T))
        ls1 = refill(800.0)
        shapes1 = [l.shape for l in ls1]
        amps1   = [l.amplitude for l in ls1]
        ls2 = refill(1600.0)
        @test ls2 === ls1
        @test [l.shape for l in ls2] == shapes1     # shapes untouched
        @test [l.amplitude for l in ls2] != amps1   # amplitudes track params
    end

    @testset "a species shared by two models keeps their instruments apart" begin
        # The static build folds the instrument into the cached shapes, so the
        # cache must belong to the model, not the species.
        species = SpeciesModel(; transitions, distribution_function = bd, lineshape_function = lsf_sta)
        broad = AbsorptionCoefficientSpectrum(; species, instrument_function = (wn, p) -> Gaussian(0.5))
        plain = AbsorptionCoefficientSpectrum(; species)
        ref   = AbsorptionCoefficientSpectrum(; transitions, distribution_function = bd, lineshape_function = lsf_dyn)
        p = ComponentVector(T = 900.0)
        broad(xs, p)
        @test plain(xs, p) ≈ ref(xs, p) rtol = 1e-13
        @test maximum(broad(xs, p)) < 0.5 * maximum(plain(xs, p))
    end
end
