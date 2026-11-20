"""
    HITRAN

Reader for HITRAN/HITEMP line lists (fixed-width "HITRAN160" `.par` format).

Covers local `.par` files (HITRAN and HITEMP share the record format), the
bundled isotopologue metadata table (`molparam`), the TIPS-2025 partition
sums (`partition_function`), and the line intensity at temperature
(`linestrength`). Not implemented: downloading, non-LTE quantum-number
decoding.
"""
module HITRAN

# `Line` and `is_usable` are not exported. `Line` would clash with
# `SpectralLines.Line`.
export read_par, read_hitemp, molecule_name, molparam, MOLECULES
export partition_function, PartitionFunction, TIPS, linestrength

include("format.jl")
include("parse.jl")
include("molparam.jl")
include("tips.jl")

end # module HITRAN
