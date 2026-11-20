using SpectraUtils
using Test
using ComponentArrays
using ForwardDiff

@testset "ChebyshevBackground" begin
    lo, hi = 300.0, 320.0
    xs = collect(range(lo, hi; length = 41))
    t(x) = (2x - (hi + lo)) / (hi - lo)     # map onto [-1, 1]
    # T0 = 1, T1 = t, T2 = 2t^2 - 1, T3 = 4t^3 - 3t
    basis = (x -> one(x), x -> t(x), x -> 2t(x)^2 - 1, x -> 4t(x)^3 - 3t(x))
    bg = ChebyshevBackground(lo, hi, 3)

    @testset "matches the Chebyshev series" begin
        c = [0.7, -0.3, 0.2, 0.05]
        expect = [sum(c[k] * basis[k](x) for k in 1:4) for x in xs]
        @test bg(xs, ComponentVector(background = c)) ≈ expect atol = 1e-13
        # Degree 0 skips the Clenshaw loop entirely.
        b0 = ChebyshevBackground(lo, hi, 0)
        @test all(≈(2.5), b0(xs, ComponentVector(background = [2.5])))
    end

    @testset "construction and a custom key" begin
        @test_throws ArgumentError ChebyshevBackground(320.0, 300.0, 2)   # hi <= lo
        @test_throws ArgumentError ChebyshevBackground(lo, hi, -1)        # negative degree
        @test_throws DimensionMismatch bg(xs, ComponentVector(background = [1.0, 2.0]))

        bk = ChebyshevBackground(lo, hi, 1; key = :continuum)
        p = ComponentVector(T = 500.0, continuum = [1.0, 0.25])
        @test bk(xs, p) ≈ ChebyshevBackground(lo, hi, 1)(xs, ComponentVector(background = [1.0, 0.25]))
    end

    @testset "works as a spectrum model's background field" begin
        tr = ElectricDipoleTransition(State(:X, 0, nothing, 0.0, 1.0),
                                      State(:X, 0, nothing, 1000.0, 3.0), 1000.0, 1.0e6)
        ν = collect(range(995.0, 1005.0; length = 201))
        bgm = ChebyshevBackground(995.0, 1005.0, 2)
        p = ComponentVector(T = 800.0, background = [0.3, 0.1, -0.05])

        mk(b) = EmissionSpectrum(;
            transitions = [tr], distribution_function = VariableDistribution(:T),
            lineshape_function = (t_, q) -> Lorentzian(0.5), background = b)

        with, without = mk(bgm)(ν, p), mk(nothing)(ν, p)
        # The difference is the background. Tolerance is set by cancellation in
        # the subtraction: the spectrum peaks near 5e8, so one ulp is ~6e-8.
        @test maximum(abs, (with .- without) .- bgm(ν, p)) < 10 * eps(maximum(without))

        ret = similar(with); mk(bgm)(ret, ν, p)
        @test ret ≈ with rtol = 1e-12
    end

    @testset "differentiable in its coefficients" begin
        f(c) = sum(bg(xs, ComponentVector(background = c)))
        g = ForwardDiff.gradient(f, [0.3, 0.1, -0.05, 0.02])
        # d/dc_k of the sum is the sum of T_k over the grid.
        @test g ≈ [sum(basis[k].(xs)) for k in 1:4] rtol = 1e-10
    end
end
