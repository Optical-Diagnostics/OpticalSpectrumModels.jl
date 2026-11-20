"""
    ChebyshevBackground(lo, hi, degree; key = :background)

A Chebyshev-polynomial background, callable as `(xs, params) -> values`, the
form every model's `background` field expects:

```julia
sp = EmissionSpectrum(; transitions, distribution_function, lineshape_function,
                        background = ChebyshevBackground(300.0, 320.0, 3))
```

The coefficients live in `params` under `key`, so they are fit parameters like
any other:

```julia
params = ComponentVector(T = 1000.0, background = [0.4, 0.05, -0.01, 0.002])
priors = ComponentVector(T = Normal(900.0, 200.0),
                         background = [Normal(0.0, 1.0) for _ in 1:4])
```

`degree + 1` coefficients are expected; [`n_coefficients`](@ref) returns that.
Degree 0 is a constant offset, degree 1 a straight line.

Give
the actual fitted range as `lo`/`hi`: points outside still evaluate, but the
conditioning holds only inside.
"""
struct ChebyshevBackground{K,T}
    lo::T
    hi::T
    degree::Int

    function ChebyshevBackground{K,T}(lo::T, hi::T, degree::Int) where {K,T}
        hi > lo || throw(ArgumentError("ChebyshevBackground: need hi > lo, got lo=$lo, hi=$hi"))
        degree >= 0 || throw(ArgumentError("ChebyshevBackground: degree must be ≥ 0, got $degree"))
        return new{K,T}(lo, hi, degree)
    end
end

function ChebyshevBackground(lo, hi, degree::Integer; key::Symbol = :background)
    T = promote_type(typeof(lo), typeof(hi), Float64)
    return ChebyshevBackground{key,T}(T(lo), T(hi), Int(degree))
end

"""
    n_coefficients(bg::ChebyshevBackground) -> Int

Number of coefficients the background expects in `params`: `degree + 1`.
"""
n_coefficients(bg::ChebyshevBackground) = bg.degree + 1

"""
    coefficient_key(::ChebyshevBackground) -> Symbol

The `params` key the coefficients are read from.
"""
coefficient_key(::ChebyshevBackground{K}) where {K} = K

# Map [lo, hi] onto Chebyshev's natural [-1, 1].
@inline _cheb_t(bg::ChebyshevBackground, x) = (2x - (bg.hi + bg.lo)) / (bg.hi - bg.lo)

# Σ cₖ Tₖ(t) by the Clenshaw recurrence.
@inline function _clenshaw(c, t)
    n = length(c)
    n == 0 && return zero(t) * zero(eltype(c))
    @inbounds begin
        n == 1 && return c[1] * one(t)
        b1 = zero(t) * zero(eltype(c))
        b2 = b1
        t2 = 2t
        for k in n:-1:2
            b1, b2 = c[k] + t2 * b1 - b2, b1
        end
        return c[1] + t * b1 - b2
    end
end

function (bg::ChebyshevBackground{K})(xs, params) where {K}
    c = params[Val(K)]
    length(c) == bg.degree + 1 || throw(DimensionMismatch(
        "ChebyshevBackground of degree $(bg.degree) expects $(bg.degree + 1) " *
        "coefficients in params[:$K], got $(length(c))"))
    return map(x -> _clenshaw(c, _cheb_t(bg, x)), xs)
end

function Base.show(io::IO, bg::ChebyshevBackground{K}) where {K}
    print(io, "ChebyshevBackground(", bg.lo, ", ", bg.hi, ", ", bg.degree, "; key=:", K, ")")
end
