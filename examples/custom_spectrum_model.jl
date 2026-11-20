using SpectraUtils
using ComponentArrays
using ForwardDiff

# A custom spectrum model built from the public API, showing what
# `LineCache`/`fill_lines!` do inside EmissionSpectrum/RadianceSpectrum.
#
# A fit calls the model with `params::ComponentVector{T}`, where `T` is
# `Float64` on one call and `ForwardDiff.Dual{...}` on the next. One
# `Vector{Line{T,...}}` cannot hold both, so the cache keeps one array per `T`.
#
# The model is minimal: emission scaled by one global "brightness" parameter,
# without a Boltzmann distribution.

struct MyEmissionSpectrum{TR}
    transitions::Vector{TR}
    line_cache::LineCache
end
MyEmissionSpectrum(transitions) = MyEmissionSpectrum(transitions, LineCache())

function (sp::MyEmissionSpectrum)(xs, params::ComponentVector{T}) where {T}
    build_line(t, p) = emission_line(t,
        (state, p) -> p.brightness,
        (t, p) -> Gaussian(p.sigma),
        p)
    instrument_function(wn, p) = nothing

    ls = fill_lines!(sp.line_cache, :emission, sp.transitions, build_line, instrument_function, params)
    return calculate_lines(xs, ls, PointSample())
end

# Synthetic transitions, so the example runs without ExoMol data
lower0 = State(:X, 0, nothing, 0.0, 1.0)
upper0 = State(:X, 0, nothing, 1000.0, 3.0)
t1     = ElectricDipoleTransition(lower0, upper0, 1000.0, 10.0)

lower1 = State(:X, 1, nothing, 0.0, 1.0)
upper1 = State(:X, 1, nothing, 1050.0, 5.0)
t2     = ElectricDipoleTransition(lower1, upper1, 1050.0, 8.0)

sp = MyEmissionSpectrum([t1, t2])
xs = collect(range(950.0, 1100.0; length=500))

params = ComponentVector(brightness=1.0, sigma=0.5)
spectrum1 = sp(xs, params)
spectrum2 = sp(xs, params)
@assert spectrum1 == spectrum2
@info "Float64 evaluation" peak=maximum(spectrum1)

f(p) = sp(xs, p)
J = ForwardDiff.jacobian(f, params)
@assert all(isfinite, J)
@info "ForwardDiff jacobian computed" size=size(J) max_abs=maximum(abs, J)

spectrum3 = sp(xs, params)
@assert spectrum3 == spectrum1
