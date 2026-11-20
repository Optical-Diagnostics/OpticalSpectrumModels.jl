"""
    State{E,V,J,EN,D}

A quantum state: `electronic_state`, `vibrational_state`, `rotational_state`
(any label type, e.g. `Symbol`/`Int`, or `nothing` if unspecified), plus optional
`energy` [cm⁻¹] and `degeneracy`. Construct with `State(electronic,
vibrational=nothing, rotational=nothing)`, which leaves energy and degeneracy
`nothing`, or with all five: `State(electronic, vibrational, rotational, energy,
degeneracy)`.

`child in parent` holds when the electronic states are equal and `parent`'s
vibrational and rotational states are each `nothing` (a wildcard) or equal to
`child`'s. `Boltzmann` tests membership in its `parent` this way;
`StateResolved` groups states by e.g. `vibrational_state`.
"""
struct State{E,V,J,EN,D}
    electronic_state::E
    vibrational_state::V
    rotational_state::J
    energy::EN
    degeneracy::D
end

State(electronic_state, vibrational_state=nothing, rotational_state=nothing) =
    State{typeof(electronic_state), typeof(vibrational_state), typeof(rotational_state), Nothing, Nothing}(electronic_state, vibrational_state, rotational_state, nothing, nothing)

Base.in(child::State, parent::State) = parent.electronic_state == child.electronic_state &&
                                       (parent.vibrational_state == child.vibrational_state || parent.vibrational_state === nothing) &&
                                       (parent.rotational_state == child.rotational_state || parent.rotational_state === nothing)

function Base.show(io::IO, s::State)
    print(io, "State(", s.electronic_state, ", ", s.vibrational_state, ", ", s.rotational_state)
    extras = String[]
    isnothing(s.energy) || push!(extras, "E=$(s.energy)")
    isnothing(s.degeneracy) || push!(extras, "g=$(s.degeneracy)")
    isempty(extras) || print(io, "; ", join(extras, ", "))
    print(io, ")")
end


"""
    upper_manifold(transitions) -> State

The electronic manifold shared by the upper states of `transitions`, as a `State`
for use as a `Boltzmann` `parent`. Throws if they span more than one manifold.

`Boltzmann(data, parent)` picks its partition function by whether `parent`
restricts the state set, so a `parent` that misdescribes the list gives the
wrong Q(T).

```julia
tl = filter(t -> electronic_state(t.upper) == :A, make_transition_list(data))
parent = upper_manifold(tl)
bd = Boltzmann(data, parent)
```
"""
function upper_manifold(transitions)
    isempty(transitions) && throw(ArgumentError("no transitions"))
    e = electronic_state(first(transitions).upper)
    for t in transitions
        electronic_state(t.upper) == e || throw(ArgumentError(
            "upper states span more than one electronic manifold ($e and " *
            "$(electronic_state(t.upper))); filter to one before taking a parent, " *
            "or build the `State` yourself if you meant to mix them"))
    end
    return State(e)
end

"""
    energy(s) -> E
    energy(s1, s2) -> E₁ - E₂

The energy of state `s` [cm⁻¹], or the difference between two states.
"""
energy(s::S) where {S<:State} = s.energy
energy(s1::S, s2::S) where {S<:State} = energy(s1) - energy(s2)

"""
    degeneracy(s) -> g

The statistical weight of state `s`.
"""
degeneracy(s::S) where {S<:State} = s.degeneracy

"""
    electronic_state(s)

The electronic label of state `s`.
"""
electronic_state(s::S) where {S<:State} = s.electronic_state

"""
    vibrational_state(s)

The vibrational label of state `s`, `nothing` where unspecified. Usable as the
`state_type` of a [`StateResolved`](@ref), to group states by vibrational level.
"""
vibrational_state(s::S) where {S<:State} = s.vibrational_state

"""
    vibronic_state(s) -> (electronic, vibrational)

The electronic and vibrational labels of state `s` as a tuple, to group states
by vibronic level in a [`StateResolved`](@ref).
"""
vibronic_state(s::S) where {S<:State} = (s.electronic_state, s.vibrational_state)

"""
    rotational_state(s)

The rotational label of state `s`, `nothing` where unspecified.
"""
rotational_state(s::S) where {S<:State} = s.rotational_state

"""
    ElectricDipoleTransition(lower, upper, wavenumber, A, pressure_shift = 0.0)

A radiative transition between a `lower` and `upper` [`State`](@ref), at
`wavenumber` [cm⁻¹] with Einstein A coefficient `A` [s⁻¹]. Construct with
`ElectricDipoleTransition(lower, upper, A)`, taking `wavenumber` from the states'
energy difference, or pass `wavenumber` directly.

The input to the spectrum models; [`make_transition_list`](@ref) builds a
`Vector` of these from a line-list source (ExoMol, or HITRAN via
[`SpectraUtils.HITRAN`](@ref)).

`pressure_shift` is the shift of the line centre per unit pressure in
cm⁻¹/atm (HITRAN's `δ_air`). The models apply it only through a
`lineshift_function` of the `SpeciesModel`, e.g. `(t, p) -> t.pressure_shift * p.P`.
"""
struct ElectricDipoleTransition{SL<:State,SU<:State,W,A,S}
    lower::SL
    upper::SU
    wavenumber::W
    A::A
    pressure_shift::S
end

ElectricDipoleTransition(lower, upper, wavenumber, A) =
    ElectricDipoleTransition(lower, upper, wavenumber, A, 0.0)
ElectricDipoleTransition(lower, upper, A) = ElectricDipoleTransition(lower, upper, energy(upper, lower), A)

function Base.show(io::IO, tr::ElectricDipoleTransition)
    print(io, "ElectricDipoleTransition(")
    show(io, tr.lower)
    print(io, " --> ")
    show(io, tr.upper)
    print(io, "; ν̃=", tr.wavenumber, ", A=", tr.A, ")")
end

"""
    absorption_cross_section(tr::ElectricDipoleTransition) -> Float64

Integrated absorption cross section in cm (equivalently cm²·cm⁻¹): the prefactor
that, times a lineshape of unit area over wavenumber (`∫φ(ν̃)dν̃ = 1`, so `φ` is
itself in cm) and a number density [cm⁻³], gives an absorption coefficient in
cm⁻¹.

From the Einstein A coefficient via the B-coefficient relations:

```math
\\int_{-\\infty}^{\\infty} \\sigma(\\tilde{\\nu})\\,d\\tilde{\\nu}
= \\frac{g_u}{g_l}\\frac{A_{ul}}{8\\pi\\,c\\,\\tilde{\\nu}_0^2}
```

with ``g_u``, ``g_l`` the upper/lower degeneracies, ``A_{ul}`` the Einstein A
coefficient [s⁻¹], ``c`` the speed of light in cm s⁻¹, and ``\\tilde{\\nu}_0`` the
transition wavenumber [cm⁻¹].
"""
absorption_cross_section(tr::T) where {T<:ElectricDipoleTransition} =
    degeneracy(tr.upper) / degeneracy(tr.lower) * tr.A / (8π * c0_cm * tr.wavenumber^2)


