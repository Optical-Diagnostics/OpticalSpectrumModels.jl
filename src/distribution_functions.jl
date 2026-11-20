"""
    DistributionFunction

Abstract supertype for state populations: callable as `f(state, params)`,
returning the population of `state` (a fraction, or a number density if the
distribution is scaled to one) at the parameters `params`. The spectrum models
call it on the upper and lower state of every transition.

Any `(state, params) -> population` function works in place of a subtype;
`DistributionFunction(f)` passes functions and subtypes through and wraps a
number in a [`ConstantDistribution`](@ref).
"""
abstract type DistributionFunction end



degeneracy(s) = error("degeneracy() not defined for $(typeof(s))")
energy(s) = error("energy() not defined for $(typeof(s))")
rotational_state(s) = error("rotational_state() not defined for $(typeof(s))")
vibrational_state(s) = error("vibrational_state() not defined for $(typeof(s))")
electronic_state(s) = error("electronic_state() not defined for $(typeof(s))")
vibronic_state(s) = error("vibronic_state() not defined for $(typeof(s))")

DistributionFunction(f::F) where {F<:DistributionFunction} = f
DistributionFunction(f::F) where {F<:Function} = f
DistributionFunction(f::F) where {F<:Real} = ConstantDistribution(f)


"""
    Boltzmann(partition_function, E0, parent[, temperature_function])
    Boltzmann(; partition_function, E0, parent, temperature_function = p -> p[Val(:T)])

Thermal population of the states under `parent` (see [`State`](@ref)):

```math
f(s) = \\frac{g_s}{Q(T)} \\exp\\left(-\\frac{hc\\,(E_s - E_0)}{k_B T}\\right)
```

and zero for states outside `parent`. `partition_function` is `Q`, callable as
`Q(T)` and computed from the same energy zero `E0` [cm⁻¹]. `T` comes from
`temperature_function(params)`, by default `params.T`; capture a different key
statically (`p -> p.T_rot`).

An unwrapped `Boltzmann` is global LTE ([`is_LTE`](@ref)), which gives
`RadianceSpectrum` its Planck fast path. With `using ExoMol`,
`Boltzmann(data, parent)` builds one from an ExoMol dataset; for HITRAN use
`HITRAN.partition_function` as `Q` with `E0 = 0`.
"""
@kwdef struct Boltzmann{PF,E,S<:State,TF} <: DistributionFunction
    partition_function::PF
    E0::E
    parent::S
    temperature_function::TF = (p) -> p[Val(:T)]
end

Boltzmann(partition_function, E0, parent::State) =
    Boltzmann(partition_function, E0, parent, (p) -> p[Val(:T)])

function (b::Boltzmann)(s, params)::eltype(params)
    s in b.parent || return zero(eltype(params))
    T = b.temperature_function(params)
    return degeneracy(s) / b.partition_function(T) * exp(-(energy(s) - b.E0) * 100 * h * c0 / kB / T)
end

# Temperature reported to RadianceSpectrum's LTE fast path. A custom
# distribution that overrides `is_LTE` to `true` must extend this too.
_temperature(b::Boltzmann, params) = b.temperature_function(params)
_temperature(df, params) = error(
    "_temperature is not defined for $(typeof(df)); an is_LTE distribution " *
    "must provide it for the Planck fast path")

function Base.show(io::IO, b::Boltzmann)
    print(io, "Boltzmann(parent=")
    show(io, b.parent)
    print(io, ", E0=", b.E0, ")")
end

"""
    is_LTE(::DistributionFunction) -> Bool

Whether a distribution function is a single global thermal population at one
temperature: `true` only for a bare, unwrapped `Boltzmann`. `RadianceSpectrum`
uses it to dispatch to the closed-form Planck source function.

A `Boltzmann` inside a `ProductDistribution`/`StateResolved` (e.g. one
rotational Boltzmann per vibrational level) is not global LTE and gives `false`.
"""
is_LTE(::DistributionFunction) = false
is_LTE(::Boltzmann) = true



"""
    StateResolved(state_type, distribution)
    StateResolved(; state_type, distribution)

A separate distribution per group of states: `state_type(s)` gives the group key
of state `s` (e.g. [`vibrational_state`](@ref)), and `distribution` maps keys to
distributions, evaluated as `distribution[key](s, params)`. States whose key is
missing have zero population. Wrap a fixed value in a
[`ConstantDistribution`](@ref). For a non-equilibrium population, such as one
rotational `Boltzmann` per vibrational level:

```julia
StateResolved(vibrational_state, Dict(0 => bd_v0, 1 => bd_v1))
```
"""
@kwdef struct StateResolved{ST,DS} <: DistributionFunction
    state_type::ST
    distribution::DS
end

(f::StateResolved{ST,DS})(s, params) where {ST,DS} = f.state_type(s) in keys(f.distribution) ? f.distribution[f.state_type(s)](s, params) : zero(eltype(params))

Base.show(io::IO, sr::StateResolved) = print(io, "StateResolved(", length(sr.distribution), " states)")

"""
    ProductDistribution(d1, d2, ...)

The product of several distributions' values for each state, e.g. a vibrational
population times a rotational one.
"""
struct ProductDistribution{T<:Tuple} <: DistributionFunction
    distributions::T
end

function (pd::ProductDistribution)(s, params)::eltype(params)
    return prod(f -> f(s, params), pd.distributions)
end

ProductDistribution(args...) = ProductDistribution(args)

function Base.show(io::IO, pd::ProductDistribution)
    print(io, "ProductDistribution(")
    for (i, d) in enumerate(pd.distributions)
        i > 1 && print(io, ", ")
        show(io, d)
    end
    print(io, ")")
end


"""
    ConstantDistribution(value)

The same population `value` for every state, independent of `params`.
"""
struct ConstantDistribution{T<:Real} <: DistributionFunction
    value::T
end
(cd::ConstantDistribution)(s, params=nothing) = cd.value
Base.show(io::IO, cd::ConstantDistribution) = print(io, "ConstantDistribution(value=", cd.value, ")")

"""
    VariableDistribution(key)

The same population for every state, read from the scalar parameter
`params[key]`, so it can be fitted. Combine with [`StateResolved`](@ref) for one
free population per group of states.
"""
struct VariableDistribution <: DistributionFunction
    key::Symbol
end
Base.show(io::IO, vd::VariableDistribution) = print(io, "VariableDistribution(key=:", vd.key, ")")


# Reads `params[sym]` through its integer index, which does not allocate.
function _index_allocation_free(params, sym)
    idx = first(getaxes(params))[sym].idx
    return params[idx]
end

(vd::VariableDistribution)(s, params)::eltype(params) = _index_allocation_free(params, vd.key)::eltype(params)


