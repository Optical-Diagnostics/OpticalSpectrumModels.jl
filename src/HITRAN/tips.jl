const _TIPS_FILE = joinpath(@__DIR__, "data", "tips2025.txt")

"""
    PartitionFunction

Total internal partition sum `Q(T)` of one isotopologue from the TIPS-2025
tables, tabulated at `T = 1, 10, 20, …, Tmax` K. Callable as `Q(T)`; values
agree with HAPI's `partitionSum`. Temperatures outside `1 ≤ T ≤ Tmax` raise an
error. Get one with [`partition_function`](@ref).

Reference: Gamache et al., "Partition sums for molecules and their
isotopologues for HITRAN2024", JQSRT 345, 109568 (2025).
doi:10.1016/j.jqsrt.2025.109568

Fields: `mol_id`, `iso`, `Tmax`, `Q` (values on the grid).
"""
struct PartitionFunction
    mol_id::Int
    iso::Int
    Tmax::Float64
    Q::Vector{Float64}
end

# Grid point i (1-based) sits at T = 1 for i = 1 and T = 10 (i - 1) after.
_tips_T(i::Int) = i == 1 ? 1.0 : 10.0 * (i - 1)

temperature_range(pf::PartitionFunction) = (1.0, pf.Tmax)

function (pf::PartitionFunction)(T)
    Q = pf.Q
    n = length(Q)
    1 <= T <= pf.Tmax || throw(DomainError(T,
        "TIPS partition sum of molecule $(pf.mol_id), isotopologue $(pf.iso) is tabulated for 1 K to $(pf.Tmax) K"))
    # First grid point at or above T, as in TIPS. `i` is at least 2 for T > 1.
    i = T <= 1.0 ? 2 : max(2, ceil(Int, T / 10) + 1)
    if i < 3 || i == n
        # Three points near the ends of the table.
        j = i < 3 ? 3 : n
        x0, x1, x2 = _tips_T(j - 2), _tips_T(j - 1), _tips_T(j)
        y0, y1, y2 = Q[j - 2], Q[j - 1], Q[j]
        return y0 * (T - x1) * (T - x2) / ((x0 - x1) * (x0 - x2)) +
               y1 * (T - x0) * (T - x2) / ((x1 - x0) * (x1 - x2)) +
               y2 * (T - x0) * (T - x1) / ((x2 - x0) * (x2 - x1))
    end
    # Four points: two below T, the first at or above, and the next.
    x0, x1, x2, x3 = _tips_T(i - 2), _tips_T(i - 1), _tips_T(i), _tips_T(i + 1)
    y0, y1, y2, y3 = Q[i - 2], Q[i - 1], Q[i], Q[i + 1]
    return y0 * (T - x1) * (T - x2) * (T - x3) / ((x0 - x1) * (x0 - x2) * (x0 - x3)) +
           y1 * (T - x0) * (T - x2) * (T - x3) / ((x1 - x0) * (x1 - x2) * (x1 - x3)) +
           y2 * (T - x0) * (T - x1) * (T - x3) / ((x2 - x0) * (x2 - x1) * (x2 - x3)) +
           y3 * (T - x0) * (T - x1) * (T - x2) / ((x3 - x0) * (x3 - x1) * (x3 - x2))
end

Base.show(io::IO, pf::PartitionFunction) =
    print(io, "PartitionFunction(molecule ", pf.mol_id, ", isotopologue ", pf.iso, ", 1-", pf.Tmax, " K)")

function _load_tips(path::AbstractString=_TIPS_FILE)
    table = Dict{Tuple{Int,Int},PartitionFunction}()
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        parts = split(s)
        mol_id = parse(Int, parts[1])
        iso = parse(Int, parts[2])
        Tmax = parse(Float64, parts[3])
        Q = [parse(Float64, p) for p in @view parts[4:end]]
        length(Q) == round(Int, Tmax / 10) + 1 || error("TIPS table for $mol_id/$iso has $(length(Q)) values, expected $(round(Int, Tmax / 10) + 1)")
        all(iszero, Q) && continue  # placeholder without data, e.g. atomic oxygen
        table[(mol_id, iso)] = PartitionFunction(mol_id, iso, Tmax, Q)
    end
    return table
end

"""
    TIPS

Partition sums of the HITRAN isotopologues, keyed by `(mol_id, iso)`, from
TIPS-2025 as distributed with HAPI 1.3 (MIT licence, HITRAN team).
Isotopologues without TIPS data, such as atomic oxygen, are absent. The
`Q(296 K)` column of [`MOLPARAM`](@ref) comes from an earlier edition and
differs by up to a factor of two for a few minor species.
"""
const TIPS = _load_tips()

"""
    partition_function(mol_id, iso) -> PartitionFunction

The TIPS-2025 partition sum `Q(T)` of a HITRAN isotopologue, callable as
`Q(T)`. Use it as the partition function of a `Boltzmann` distribution over
HITRAN transitions:

```julia
Q = HITRAN.partition_function(2, 1)              # CO2, 626
bd = Boltzmann(Q, 0.0, State(:X))                 # energies count from the ground state
```

`Q(T)` includes the nuclear spin degeneracy, as do the HITRAN `g'`, `g''`
values, so populations from the two are consistent.
"""
function partition_function(mol_id::Integer, iso::Integer)
    haskey(TIPS, (Int(mol_id), Int(iso))) ||
        throw(KeyError("no TIPS partition sum for molecule $mol_id, isotopologue $iso"))
    return TIPS[(Int(mol_id), Int(iso))]
end

# Second radiation constant h c / k in cm K, CODATA 2018 exact values.
const C2_CM_K = 6.62607015e-34 * 2.99792458e10 / 1.380649e-23

"""
    linestrength(S_ref, El, wav, Q, T; Tref = 296.0)
    linestrength(line::Line, T; Tref = 296.0)

HITRAN line intensity at temperature `T` from its value `S_ref` at `Tref`:

    S(T) = S_ref · Q(Tref)/Q(T) · exp(-c₂ El / T) / exp(-c₂ El / Tref)
                 · (1 - exp(-c₂ ν / T)) / (1 - exp(-c₂ ν / Tref))

with `El` the lower-state energy and `wav` the wavenumber in cm⁻¹, `Q` the
partition sum callable as `Q(T)`, and `c₂ = hc/k`. The `Line` method takes
`intensity`, `El`, `wav` from the line and `Q` from [`partition_function`](@ref).
Like HITRAN's `S`, the result includes the natural isotopic abundance.
"""
function linestrength(S_ref, El, wav, Q, T; Tref = 296.0)
    boltzmann = exp(-C2_CM_K * El * (1 / T - 1 / Tref))
    emission = expm1(-C2_CM_K * wav / T) / expm1(-C2_CM_K * wav / Tref)
    return S_ref * (Q(Tref) / Q(T)) * boltzmann * emission
end

linestrength(l::Line, T; Tref = 296.0) =
    linestrength(l.intensity, l.El, l.wav, partition_function(l.mol_id, l.iso), T; Tref)
