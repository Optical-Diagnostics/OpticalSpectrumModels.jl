const _MOLPARAM_FILE = joinpath(@__DIR__, "data", "molparam.txt")

const MolParamEntry = @NamedTuple{isoname::String, abundance::Float64, Q296K::Float64, gj::Int, molar_mass::Float64}

function _load_molparam(path::AbstractString=_MOLPARAM_FILE)
    table = Dict{Tuple{Int,Int},MolParamEntry}()
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        parts = split(s)
        length(parts) == 7 || continue
        id = tryparse(Int, parts[1])
        id === nothing && continue      # header row ("id iso isoname ...")
        iso = parse(Int, parts[2])
        table[(id, iso)] = (
            isoname = String(parts[3]),
            abundance = parse(Float64, parts[4]),
            Q296K = parse(Float64, parts[5]),
            gj = parse(Int, parts[6]),
            molar_mass = parse(Float64, parts[7]),
        )
    end
    return table
end

"""
    MOLPARAM

Isotopologue metadata keyed by `(mol_id, iso)`: natural abundance, `Q(296 K)`,
nuclear-spin degeneracy `gj`, molar mass. From HITRAN's public `molparam.txt`
(http://hitran.org/media/molparam.txt).
"""
const MOLPARAM = _load_molparam()

"""
    molparam(mol_id, iso) -> (; isoname, abundance, Q296K, gj, molar_mass)

Look up isotopologue metadata for a given HITRAN molecule/isotopologue pair.
"""
molparam(mol_id::Integer, iso::Integer) = MOLPARAM[(Int(mol_id), Int(iso))]
