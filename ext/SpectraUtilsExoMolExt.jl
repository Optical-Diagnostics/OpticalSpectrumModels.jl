module SpectraUtilsExoMolExt

using ExoMol
import SpectraUtils: State, ElectricDipoleTransition, Boltzmann, h, c0, kB, u_kg, molecular_mass_kg

import SpectraUtils: make_transition_list
import SpectraUtils
using DataInterpolations


struct ExoMolState{T}
    state::T
end

# ExoMol datasets use different column names per molecule.
function State(exo::ExoMolState)
    s = exo.state
    elec = hasproperty(s, :ElecState) ? Symbol(s.ElecState) : :X
    v    = hasproperty(s, :v)               ? Int8(s.v) :
           hasproperty(s, Symbol("dcs:v"))  ? Int8(s.var"dcs:v") : nothing
    # J is half-integer for some molecules (OH's X²Π/A²Σ⁺): Int32 when whole,
    # Float64 otherwise.
    j    = hasproperty(s, :J) ? (isinteger(s.J) ? Int32(s.J) : Float64(s.J)) : nothing
    g    = hasproperty(s, :gtot)  ? Int32(s.gtot) :
           hasproperty(s, :g_tot) ? Int32(s.g_tot) : nothing
    State(elec, v, j, s.E, g)
end

"""
    molecular_mass_kg(iso::ExoMol.Isotopologue) -> m [kg]

Isotopologue mass from the dataset's own `.def` file, in kilograms, so that
`sigma_doppler(t.wavenumber, T, iso)` works directly.
"""
molecular_mass_kg(iso::ExoMol.Isotopologue) = ExoMol.molar_mass(iso) * u_kg

"""
    make_transition_list(data::ExoMol.Isotopologue) -> Vector{ElectricDipoleTransition}

One transition per line of the dataset, with upper and lower [`State`](@ref)s
built from its states table: electronic state, `v`, `J`, energy and total
degeneracy, under whichever column names the dataset uses.
"""
function make_transition_list(exo_data::I) where {I <: ExoMol.Isotopologue}
    # State-ID field name, detected once: :ID (WCCRMT) or :StateID (EBJT, …).
    sample = first(exo_data.states)
    id_sym = hasproperty(sample, :ID) ? :ID : :StateID

    # The transitions index states by ID, so the states table must be sorted
    # by it.
    function checked_state(id)
        s = exo_data.states[id]
        id == getproperty(s, id_sym) || throw(ArgumentError(
            "state $id of the states table has ID $(getproperty(s, id_sym)); " *
            "the table must be sorted by ID with no gaps"))
        return ExoMolState(s)
    end
    function make_transition(t)
        lower = checked_state(t.lower_id)
        upper = checked_state(t.upper_id)
        return ElectricDipoleTransition(State(lower), State(upper), t.wavenumber, t.A)
    end
    return map(make_transition, exo_data.transitions)
end


# Default temperature grid of the state-sum spline, covering plasma temperatures.
const _DEFAULT_TS = 50:10:10000

"""
    StateSum

Partition function Q(T) = Σ g exp(−hc (E − E0) / kT) over a set of states: a cubic
spline on a temperature grid for speed, and the sum itself outside the grid.
"""
struct StateSum{S}
    spline::S
    g::Vector{Float64}       # degeneracies
    e::Vector{Float64}       # energies above E0 times hc/k, in K
    Tmin::Float64
    Tmax::Float64
end

_state_sum(g, e, T) = sum(i -> g[i] * exp(-e[i] / T), eachindex(g, e))

(q::StateSum)(T) = q.Tmin <= T <= q.Tmax ? q.spline(T) : _state_sum(q.g, q.e, T)

function partition_function(states, parent; Ts = _DEFAULT_TS)
    included_states = filter(s -> s in parent, states)
    E0 = minimum(s -> s.energy, included_states)
    g = Float64[s.degeneracy for s in included_states]
    e = Float64[(s.energy - E0) * 100 * h * c0 / kB for s in included_states]
    spline = CubicSpline([_state_sum(g, e, T) for T in Ts], Ts)
    return StateSum(spline, g, e, Float64(first(Ts)), Float64(last(Ts))), E0
end

"""
    Boltzmann(data::ExoMol.Isotopologue, parent::State; Ts = nothing,
              temperature_function = p -> p[Val(:T)])

Thermal-equilibrium distribution over the states of `data` falling under
`parent`. `temperature_function` is as for [`Boltzmann`](@ref).

If every loaded state is in `parent` and the dataset has a `.pf` file, Q(T)
comes from that file. It is the total over all states, relative to E = 0, and
more accurate than a sum over an intensity-truncated states file.

Otherwise `parent` selects a subset, such as one electronic state of OH
`MYTHOS` or one vibrational level for a `StateResolved` population. Q(T) is
then the sum over that subset. `Ts` applies only here: it is the grid of a
spline that speeds up the sum, by default 50–10 000 K in 10 K steps. Outside
the grid the sum is evaluated directly.
"""
function SpectraUtils.Boltzmann(data::ExoMol.Isotopologue, parent::State; Ts = nothing,
                                temperature_function = (p) -> p[Val(:T)])
    states = State.(ExoMolState.(data.states))
    if data.partition_function !== nothing && all(s -> s in parent, states)
        Ts === nothing || @warn(
            "Boltzmann: `Ts` is ignored, since Q(T) comes from the dataset's " *
            ".pf file. `Ts` applies only when `parent` restricts the state set.",
            Ts, parent, maxlog = 1)
        # E0 = 0 is the energy zero of the `.pf` file. The state-sum path
        # measures from the lowest included state.
        return Boltzmann(data.partition_function, zero(first(states).energy), parent, temperature_function)
    end
    return Boltzmann(partition_function(states, parent; Ts = something(Ts, _DEFAULT_TS))..., parent, temperature_function)
end

end
