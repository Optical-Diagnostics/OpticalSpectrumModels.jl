"""
    Line

One parsed HITRAN/HITEMP transition record ("HITRAN160" `.par` format). Field
names and units follow the HITRAN convention:

- `mol_id`, `iso`: HITRAN molecule and isotopologue identifiers
- `wav`: vacuum wavenumber [cm⁻¹]
- `intensity`: line intensity at 296 K [cm⁻¹/(molecule·cm⁻²)]
- `A`: Einstein A coefficient [s⁻¹]
- `airbrd`, `selbrd`: air-/self-broadened HWHM at 296 K [cm⁻¹·atm⁻¹]
- `El`: lower-state energy [cm⁻¹]
- `Tdpair`: temperature-dependence exponent for `airbrd`
- `Pshft`: air-pressure-induced line shift at 296 K [cm⁻¹·atm⁻¹]
- `globu`, `globl`, `locu`, `locl`: raw global/local quantum-number strings
  (upper/lower); not decoded into `v`/`J`/electronic state
- `ierr`, `iref`: uncertainty and reference index codes (raw strings)
- `lmix`: line-mixing flag
- `gp`, `gpp`: upper-/lower-state statistical weights (degeneracies)

Blank numeric fields parse to `NaN`: not every molecule populates every column.
"""
struct Line
    mol_id::Int
    iso::Int
    wav::Float64
    intensity::Float64
    A::Float64
    airbrd::Float64
    selbrd::Float64
    El::Float64
    Tdpair::Float64
    Pshft::Float64
    globu::String
    globl::String
    locu::String
    locl::String
    ierr::String
    iref::String
    lmix::Char
    gp::Float64
    gpp::Float64
end

"""
    is_usable(line::Line) -> Bool

Whether `line` carries what the spectrum models need: a positive Einstein `A`,
positive statistical weights `gp`/`gpp`, and a non-negative lower-state energy
`El`. Unassigned lines leave these blank, zero or negative, and
`make_transition_list` rejects them. Drop them first:

```julia
lines = filter(HITRAN.is_usable, HITRAN.read_par("H2O.par"))
```
"""
is_usable(l::Line) = isfinite(l.A) && l.A > 0 && l.gp > 0 && l.gpp > 0 && l.El >= 0

Base.show(io::IO, l::Line) =
    print(io, "Line(mol_id=", l.mol_id, ", iso=", l.iso, ", wav=", l.wav, ", A=", l.A, ")")

@inline function _parsefloat(s::AbstractString)
    t = strip(s)
    isempty(t) ? NaN : parse(Float64, t)
end

@inline function _field(rec, name::Symbol)
    for (n, a, b) in HITRAN160_OFFSETS
        n === name && return String(rec[a:b])
    end
    error("Unknown HITRAN160 field $name")
end

# CO2 isotopologues 10/11/12 are encoded as "0"/"A"/"B" in the single-character
# `iso` field instead of a plain digit.
const _CO2_ISO_QUIRK = Dict("0" => 10, "A" => 11, "B" => 12)

@inline function _parse_iso(mol_id::Int, s::AbstractString)
    if mol_id == 2
        remapped = get(_CO2_ISO_QUIRK, s, nothing)
        remapped === nothing || return remapped
    end
    return parse(Int, s)
end

function _parse_line(rec, mol_id::Int, wav::Float64)
    iso = _parse_iso(mol_id, strip(_field(rec, :iso)))
    Line(
        mol_id, iso, wav,
        _parsefloat(_field(rec, :int)),
        _parsefloat(_field(rec, :A)),
        _parsefloat(_field(rec, :airbrd)),
        _parsefloat(_field(rec, :selbrd)),
        _parsefloat(_field(rec, :El)),
        _parsefloat(_field(rec, :Tdpair)),
        _parsefloat(_field(rec, :Pshft)),
        _field(rec, :globu),
        _field(rec, :globl),
        _field(rec, :locu),
        _field(rec, :locl),
        _field(rec, :ierr),
        _field(rec, :iref),
        only(_field(rec, :lmix)),
        _parsefloat(_field(rec, :gp)),
        _parsefloat(_field(rec, :gpp)),
    )
end

function _terminator_width(bytes::Vector{UInt8})
    for i in eachindex(bytes)
        if bytes[i] == UInt8('\n')
            return (i >= 2 && bytes[i - 1] == UInt8('\r')) ? 2 : 1
        end
    end
    return 1  # no line return found (single-record file)
end

"""
    read_par(file; wn_min=-Inf, wn_max=Inf) -> Vector{Line}

Read a HITRAN or HITEMP `.par` file (both use the fixed-width "HITRAN160"
format) into a `Vector{Line}`.

Only records within `[wn_min, wn_max]` are kept, which bounds memory and time
when loading a narrow window out of a large HITEMP file.
"""
function read_par(file::AbstractString; wn_min::Real=-Inf, wn_max::Real=Inf)
    bytes = read(file)
    n = length(bytes)
    n < HITRAN160_WIDTH && return Line[]

    term = _terminator_width(bytes)
    record_width = HITRAN160_WIDTH + term

    lines = Line[]
    sizehint!(lines, cld(n, record_width))

    pos = 1
    while pos + HITRAN160_WIDTH - 1 <= n
        rec = @view bytes[pos:(pos + HITRAN160_WIDTH - 1)]
        wav = _parsefloat(_field(rec, :wav))
        if wn_min <= wav <= wn_max
            # Per record: a hitran.org export can mix molecules in one file.
            mol_id = parse(Int, strip(_field(rec, :id)))
            push!(lines, _parse_line(rec, mol_id, wav))
        end
        pos += record_width
    end

    return lines
end

# HITEMP chunk-filename convention, e.g. "02_00000-00500_HITEMP2010.par":
# molecule id "02", wavenumber band 0-500 cm⁻¹. Non-matching names return `true`
# (don't skip).
function _file_overlaps_range(path::AbstractString, wn_min::Real, wn_max::Real)
    (isinf(wn_min) && isinf(wn_max)) && return true
    m = match(r"_(\d+)-(\d+)_", basename(path))
    m === nothing && return true
    lo = parse(Float64, m.captures[1])
    hi = parse(Float64, m.captures[2])
    return lo <= wn_max && hi >= wn_min
end

"""
    read_hitemp(files; wn_min=-Inf, wn_max=Inf) -> Vector{Line}
    read_hitemp(directory; wn_min=-Inf, wn_max=Inf) -> Vector{Line}

Read and concatenate HITEMP `.par` files for one molecule (an explicit list, or
every `.par` in `directory`), such as the wavenumber-banded chunks HITEMP ships
large molecules (CO2, H2O) in. Files are read in parallel across threads.

Only records within `[wn_min, wn_max]` are kept, as in [`read_par`](@ref).
"""
function read_hitemp(files::AbstractVector{<:AbstractString}; wn_min::Real=-Inf, wn_max::Real=Inf)
    relevant = filter(f -> _file_overlaps_range(f, wn_min, wn_max), files)
    chunks = Vector{Vector{Line}}(undef, length(relevant))
    @sync for (i, f) in enumerate(relevant)
        Threads.@spawn chunks[i] = read_par(f; wn_min, wn_max)
    end
    return reduce(vcat, chunks; init=Line[])
end

read_hitemp(directory::AbstractString; wn_min::Real=-Inf, wn_max::Real=Inf) =
    read_hitemp(
        filter(f -> endswith(f, ".par"), joinpath.(directory, readdir(directory)));
        wn_min, wn_max,
    )
