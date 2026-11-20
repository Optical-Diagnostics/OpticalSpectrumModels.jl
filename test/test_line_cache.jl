using SpectraUtils
using Test
using ComponentArrays
using ForwardDiff
import SpectraUtils: convolve, LineShape

include("synthetic.jl")

# A line container other than `Line`, shaped like the internal `_DualLine`.
# `fill_lines!` must build its buffer from the container type, not assume `Line`.
struct _TestDualLine{T, LS<:LineShape{T}}
    position::T
    amp1::T
    amp2::T
    shape::LS
end
function _TestDualLine(position, amp1, amp2, shape::LS) where {T, LS<:LineShape{T}}
    num_type = promote_type(typeof(position), typeof(amp1), typeof(amp2), T)
    line_type = Base.typename(LS).wrapper{num_type}
    return _TestDualLine{num_type, line_type}(position, amp1, amp2, line_type(shape))
end
convolve(dl::_TestDualLine, ls::LineShape) = _TestDualLine(dl.position, dl.amp1, dl.amp2, convolve(dl.shape, ls))

build_emission(t, p) = emission_line(t, VariableDistribution(:f), (t, p) -> Lorentzian(0.5), p)
no_instrument(wn, p) = nothing

@testset "LineCache / fill_lines!" begin
    params = ComponentVector(f = 0.7)

    @testset "reuse and resize, independently per slot" begin
        cache = LineCache()
        build_absorption(t, p) = Line(one(eltype(p)) * t.wavenumber, one(eltype(p)) * 2 * p[:f], Gaussian(0.2))

        ls_e = fill_lines!(cache, :emission, transitions, build_emission, no_instrument, params)
        ls_a = fill_lines!(cache, :absorption, transitions, build_absorption, no_instrument, params)
        @test [l.amplitude for l in ls_e] ≈ 0.7 .* [t1.A, t2.A]
        @test all(l -> l.shape isa Gaussian && l.amplitude ≈ 2 * 0.7, ls_a)

        @test fill_lines!(cache, :emission, transitions, build_emission, no_instrument, params) === ls_e
        @test fill_lines!(cache, :absorption, transitions, build_absorption, no_instrument, params) === ls_a
        # Same array, resized in place.
        ls_short = fill_lines!(cache, :emission, [t1], build_emission, no_instrument, params)
        @test ls_short === ls_e && length(ls_short) == 1
    end

    @testset "a non-Line container across a number-type change" begin
        # The buffer takes the container type of the sample line, as
        # `RadianceSpectrum`'s `_DualLine` slot needs.
        cache = LineCache()
        build_dual(t, p) = _TestDualLine(one(eltype(p)) * t.wavenumber, one(eltype(p)) * p[:f], one(eltype(p)) * 2p[:f], Lorentzian(0.5))
        params_d = ComponentVector(f = ForwardDiff.Dual(0.7, 1.0))

        ls_f = fill_lines!(cache, :dual, transitions, build_dual, no_instrument, params)
        @test all(l -> l isa _TestDualLine{Float64}, ls_f)
        ls_d = fill_lines!(cache, :dual, transitions, build_dual, no_instrument, params_d)
        @test all(l -> l isa _TestDualLine{<:ForwardDiff.Dual}, ls_d)
        # The Float64 slot survives the Dual call.
        @test fill_lines!(cache, :dual, transitions, build_dual, no_instrument, params) === ls_f
    end

    @testset "AD through a cached call" begin
        sp = EmissionSpectrum(; transitions, distribution_function = VariableDistribution(:f),
                                lineshape_function = (t, p) -> Lorentzian(0.5))
        xs = collect(range(950.0, 1100.0; length = 50))
        p1 = ComponentVector(f = 1.0)
        sp(xs, p1)   # fill the Float64 slot first
        J = ForwardDiff.jacobian(p -> sp(xs, p), p1)
        # The spectrum is linear in f, so ∂/∂f at f = 1 is the spectrum itself.
        @test J[:, 1] ≈ sp(xs, p1) rtol = 1e-12
    end
end
