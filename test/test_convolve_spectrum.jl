using SpectraUtils
using Test

@testset "convolve_spectrum" begin
    @testset "Gaussian instrument matches the analytic convolution" begin
        # Gaussian ⊛ Gaussian is Gaussian(√(0.2² + 0.5²)).
        natural = Gaussian(0.2)
        model(xs, params) = natural.(xs)
        instrument = Gaussian(0.5)
        xs = range(-10.0, 10.0; step = 0.01)
        result = SpectraUtils.convolve_spectrum(model, xs, nothing, instrument)
        @test result ≈ SpectraUtils.convolve(natural, instrument).(xs) atol = 1e-4
    end

    @testset "Triangular slit reproduces the RADIS convention (zero at |x| >= fwhm)" begin
        # A narrow Gaussian convolved with Triangular(w) is close to Triangular(w).
        natural = Gaussian(0.01)
        model(xs, params) = natural.(xs)
        xs = range(-3.0, 3.0; step = 0.001)
        result = SpectraUtils.convolve_spectrum(model, xs, nothing, Triangular(1.0))
        # The finite Gaussian rounds Triangular's corners.
        @test maximum(abs.(result .- Triangular(1.0).(xs))) < 1e-2
    end

    @testset "a slightly non-uniform grid stays registered to xs" begin
        # Spacing grows by 0.8% across the grid, inside the 1% tolerance.
        n = 1001
        steps = [0.01 * (1 + 0.008 * (i - 1) / (n - 2)) for i in 1:n-1]
        xs = cumsum(vcat(0.0, steps))
        model(x, params) = Gaussian(0.05).(x .- xs[900])
        result = SpectraUtils.convolve_spectrum(model, xs, nothing, Gaussian(0.02))
        @test argmax(result) == 900
    end

    @testset "rejects grids it cannot convolve on" begin
        model(xs, params) = Gaussian(0.2).(xs)
        xs = collect(range(-5.0, 5.0; step = 0.01))
        @test_throws ArgumentError SpectraUtils.convolve_spectrum(model, reverse(xs), nothing, Triangular(0.5))
        xn = vcat(-5.0:0.01:0.0, 0.05:0.05:5.0)
        @test_throws ArgumentError SpectraUtils.convolve_spectrum(model, xn, nothing, Triangular(0.5))
    end
end
