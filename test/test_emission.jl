using SpectraUtils
using Test
using ComponentArrays

const tr = ElectricDipoleTransition(State(:X, 0, 1), State(:A, 1, 1), 1000.0, 1e6)
const base = (
    transitions           = [tr],
    distribution_function = VariableDistribution(:T),
    lineshape_function    = (t, p) -> Lorentzian(0.5),
)

@testset "BinAverage matches PointSample on a resolving grid" begin
    xs = collect(range(995.0, 1005.0; length = 1001))
    params = ComponentVector(T = 500.0)
    ref = EmissionSpectrum(; base...)(xs, params)
    binned = EmissionSpectrum(; species = SpeciesModel(; base..., method = BinAverage(bin_width = xs[2] - xs[1])))(xs, params)
    # A bin average equals the point value up to Δx²·φ''/24, ~3e-5 at the peak.
    @test maximum(abs, binned .- ref) / maximum(ref) < 1e-4
end

@testset "number density, path length, wavenumber mapping and background" begin
    xs = collect(range(995.0, 1005.0; length = 501))
    params = ComponentVector(T = 500.0, N = 2.5, L = 3.0)
    plain = EmissionSpectrum(; base...)
    y0 = plain(xs, params)

    scaled = EmissionSpectrum(; base...,
        number_density_function = (p) -> p[:N],
        path_length_function    = (p) -> p[:L])
    @test scaled(xs, params) ≈ y0 .* (2.5 * 3.0) rtol = 1e-12

    shifted = EmissionSpectrum(; base..., wavenumber_function = (x, p) -> x .+ 0.3,
                               background = 7.0)
    ya = shifted(xs, params)
    @test ya ≈ plain(xs .+ 0.3, params) .+ 7.0 rtol = 1e-12
    rb = similar(ya); shifted(rb, xs, params)
    @test rb ≈ ya rtol = 1e-12
end
