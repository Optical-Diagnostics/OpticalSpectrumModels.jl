using SpectraUtils
using Test
using ComponentArrays

include("synthetic.jl")

const base = (
    transitions           = transitions,
    distribution_function = bd,
    lineshape_function    = (t, p) -> Lorentzian(0.5),
)
const nl = (number_density_function = (p) -> p[:N],
            path_length_function    = (p) -> p[:L])

@testset "Absorption-side models" begin
    xs     = collect(range(950.0, 1100.0; length = 400))
    params = ComponentVector(T = 800.0, N = 1e15, L = 10.0)
    T0, N, L = 800.0, 1e15, 10.0

    # α(ν) by hand, for a given line shape.
    function manual(shape)
        α = zeros(length(xs))
        for t in transitions
            f_l = t.lower.degeneracy / Qfun(T0) * exp(-c2 * t.lower.energy / T0)
            f_u = t.upper.degeneracy / Qfun(T0) * exp(-c2 * t.upper.energy / T0)
            σ0  = t.upper.degeneracy / t.lower.degeneracy * t.A / (8π * c0_cm * t.wavenumber^2)
            amp = N * σ0 * (f_l - f_u * t.lower.degeneracy / t.upper.degeneracy)
            α .+= amp .* shape.(xs .- t.wavenumber)
        end
        return α
    end
    manual_α = manual(Lorentzian(0.5))

    @testset "AbsorptionCoefficientSpectrum ≡ manual α" begin
        sp = AbsorptionCoefficientSpectrum(; base..., number_density_function = (p) -> p[:N])
        @test sp(xs, params) ≈ manual_α rtol = 1e-12
        # Default N = 1: the effective cross-section spectrum.
        @test AbsorptionCoefficientSpectrum(; base...)(xs, params) ≈ manual_α ./ N rtol = 1e-12

        # Linear model: per-line instrumental convolution is allowed and exact.
        # Lorentzian(0.5) ⊛ Gaussian(0.3) is Voigt(0.3, 0.5).
        sp_inst = AbsorptionCoefficientSpectrum(; base...,
            number_density_function = (p) -> p[:N],
            instrument_function = (wn, p) -> Gaussian(0.3))
        @test sp_inst(xs, params) ≈ manual(Voigt(0.3, 0.5)) rtol = 1e-12
    end

    @testset "AbsorbanceSpectrum = α·L (Napierian); TransmittanceSpectrum = exp(−αL)" begin
        A  = AbsorbanceSpectrum(; base..., nl...)(xs, params)
        Tr = TransmittanceSpectrum(; base..., nl...)(xs, params)
        @test A ≈ manual_α .* L rtol = 1e-12
        @test Tr ≈ exp.(-manual_α .* L) rtol = 1e-12

        # Kirchhoff consistency: LTE radiance with no incident light is
        # planck_radiance·(1 − T) on the same α.
        rad = RadianceSpectrum(; base..., nl...)(xs, params)
        @test rad ≈ planck_radiance.(xs, T0) .* (1 .- Tr) rtol = 1e-12
    end

    @testset "nonlinear transmittance rejects per-line instrument shapes" begin
        sp = TransmittanceSpectrum(; base..., nl...,
            instrument_function = (wn, p) -> Gaussian(0.3))
        @test_throws ArgumentError sp(xs, params)
    end

    @testset "in-place forms match" begin
        for sp in (AbsorptionCoefficientSpectrum(; base..., number_density_function = (p) -> p[:N]),
                   AbsorbanceSpectrum(; base..., nl...),
                   TransmittanceSpectrum(; base..., nl...))
            ref = sp(xs, params)
            ret = similar(ref)
            sp(ret, xs, params)
            sp(ret, xs, params)   # buffer reuse must not corrupt
            @test ret == ref
        end
    end
end
