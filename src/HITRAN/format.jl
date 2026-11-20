# The "HITRAN160" fixed-width `.par` record format, used by both HITRAN and
# HITEMP. Reference: Rothman et al., "The HITRAN 2004 molecular spectroscopic
# database", JQSRT 96 (2005).
const HITRAN160_COLUMNS = (
    (:id,     2),
    (:iso,    1),
    (:wav,   12),
    (:int,   10),
    (:A,     10),
    (:airbrd, 5),
    (:selbrd, 5),
    (:El,    10),
    (:Tdpair, 4),
    (:Pshft,  8),
    (:globu, 15),
    (:globl, 15),
    (:locu,  15),
    (:locl,  15),
    (:ierr,   6),
    (:iref,  12),
    (:lmix,   1),
    (:gp,     7),
    (:gpp,    7),
)

function _column_offsets(columns)
    offsets = Vector{Tuple{Symbol,Int,Int}}(undef, length(columns))
    pos = 1
    for (i, (name, width)) in enumerate(columns)
        offsets[i] = (name, pos, pos + width - 1)
        pos += width
    end
    return offsets
end

const HITRAN160_OFFSETS = _column_offsets(HITRAN160_COLUMNS)
const HITRAN160_WIDTH = sum(w for (_, w) in HITRAN160_COLUMNS)

# HITRAN molecule identifiers (61 entries), Rothman et al. 2005/2020 conventions.
const MOLECULES = Dict{Int,String}(
    1 => "H2O", 2 => "CO2", 3 => "O3", 4 => "N2O", 5 => "CO", 6 => "CH4",
    7 => "O2", 8 => "NO", 9 => "SO2", 10 => "NO2", 11 => "NH3", 12 => "HNO3",
    13 => "OH", 14 => "HF", 15 => "HCl", 16 => "HBr", 17 => "HI", 18 => "ClO",
    19 => "OCS", 20 => "H2CO", 21 => "HOCl", 22 => "N2", 23 => "HCN",
    24 => "CH3Cl", 25 => "H2O2", 26 => "C2H2", 27 => "C2H6", 28 => "PH3",
    29 => "COF2", 30 => "SF6", 31 => "H2S", 32 => "HCOOH", 33 => "HO2",
    34 => "O", 35 => "ClONO2", 36 => "NO+", 37 => "HOBr", 38 => "C2H4",
    39 => "CH3OH", 40 => "CH3Br", 41 => "CH3CN", 42 => "CF4", 43 => "C4H2",
    44 => "HC3N", 45 => "H2", 46 => "CS", 47 => "SO3", 48 => "C2N2",
    49 => "COCl2", 50 => "SO", 51 => "CH3F", 52 => "GeH4", 53 => "CS2",
    54 => "CH3I", 55 => "NF3", 56 => "H3+", 57 => "CH3", 58 => "S2",
    59 => "COFCl", 60 => "HONO", 61 => "ClNO2",
)

"""
    molecule_name(mol_id) -> String

Look up the HITRAN molecule name for a numeric molecule identifier
(the `id` field of a HITRAN160 record, e.g. `5` for `"CO"`).
"""
molecule_name(mol_id::Integer) = MOLECULES[Int(mol_id)]
