import PreallocationTools
using ComponentArrays: ComponentVector

"""
    LineCache()

Cache of per-transition line vectors, so a spectrum model can rebuild its lines
on every call (fitting, AD) without reallocating.

One `LineCache` covers any number of slots, requested by name in
[`fill_lines!`](@ref). `EmissionSpectrum`/`RadianceSpectrum` use it internally;
it is exported for custom models, see `examples/custom_spectrum_model.jl`.

A `LineCache` is not thread-safe, so use one model instance per thread. One
buffer is kept per numeric type seen and never freed.
"""
struct LineCache
    slots::Dict{Symbol, PreallocationTools.GeneralLazyBufferCache}
end
LineCache() = LineCache(Dict{Symbol, PreallocationTools.GeneralLazyBufferCache}())

Base.show(io::IO, c::LineCache) = print(io, "LineCache(", join(sort(collect(keys(c.slots))), ", "), ")")

# Grid-sized scratch vector from a cache slot, keyed by the type of the sample
# scalar `x1`.
function _scratch_vector!(cache::LineCache, name::Symbol, x1, n::Int)
    slot = get!(() -> PreallocationTools.GeneralLazyBufferCache(v -> typeof(v)[]),
                cache.slots, name)
    buf = PreallocationTools.get_tmp(slot, x1)::Vector{typeof(x1)}
    resize!(buf, n)
    return buf
end

"""
    update_or_fill_lines!(cache, name, amps_name, transitions, build_line,
                          build_amplitude, instrument_function, params) -> Vector

[`fill_lines!`](@ref) for `params`-independent shapes (see `StaticLineshapes`):
when the slot already holds a line vector of the right concrete type and length,
only the amplitudes are recomputed (`build_amplitude(t, params)`) and written
with `SpectralLines.update_lines!`, keeping the shapes and instrument convolutions of
the first build. Any type change (a `ForwardDiff` pass, a different chunk size)
keys a fresh buffer and rebuilds, as in `fill_lines!`.

Staleness is checked on the length and on the first line, whose position and
convolved shape must equal the freshly built ones; otherwise everything is
rebuilt. Do not mutate `transitions` in place beyond that.
`amps_name` is a second slot holding the amplitude scratch vector.
"""
function update_or_fill_lines!(cache::LineCache, name::Symbol, amps_name::Symbol,
                               transitions, build_line, build_amplitude,
                               instrument_function, params::ComponentVector{T}) where {T}
    isempty(transitions) && throw(ArgumentError(
        "update_or_fill_lines!: `transitions` is empty: the spectrum model needs " *
        "at least one transition (check the filter that built the list)"))
    l1 = convolve(build_line(transitions[1], params),
                  instrument_function(transitions[1].wavenumber, params))
    slot = get!(() -> PreallocationTools.GeneralLazyBufferCache(l -> typeof(l)[]),
                cache.slots, name)
    ls = PreallocationTools.get_tmp(slot, l1)::Vector{typeof(l1)}
    # `l1` is rebuilt on every call, so it shows whether the cached shapes are
    # still valid. `===` compares bits, ForwardDiff partials included.
    if length(ls) == length(transitions) &&
       ls[1].position === l1.position && ls[1].shape === l1.shape
        amps = _scratch_vector!(cache, amps_name, l1.amplitude, length(transitions))
        amps[1] = l1.amplitude
        @inbounds for i in 2:length(transitions)
            amps[i] = build_amplitude(transitions[i], params)
        end
        return update_lines!(ls; amplitudes = amps)
    end
    resize!(ls, length(transitions))
    ls[1] = l1
    @inbounds for i in 2:length(transitions)
        t = transitions[i]
        ls[i] = convolve(build_line(t, params), instrument_function(t.wavenumber, params))
    end
    return ls
end

"""
    fill_lines!(cache::LineCache, name::Symbol, transitions, build_line,
                instrument_function, params::ComponentVector) -> Vector

Build a vector of per-transition lines, reusing a cached buffer of the right type
and length. Each line is `build_line(t, params)` (numeric type from
`eltype(params)`) convolved with `instrument_function(t.wavenumber, params)`.

`name` selects the slot; different names are independent. A `Float64` call and
a ForwardDiff pass on the same slot do not collide.
"""
function fill_lines!(cache::LineCache, name::Symbol, transitions, build_line,
                      instrument_function, params::ComponentVector{T}) where {T}
    isempty(transitions) && throw(ArgumentError(
        "fill_lines!: `transitions` is empty: the spectrum model needs at " *
        "least one transition (check the filter that built the list)"))
    l1 = convolve(build_line(transitions[1], params),
                  instrument_function(transitions[1].wavenumber, params))

    # The buffer is keyed by the type of the first line.
    slot = get!(() -> PreallocationTools.GeneralLazyBufferCache(l -> typeof(l)[]),
                cache.slots, name)
    ls = PreallocationTools.get_tmp(slot, l1)::Vector{typeof(l1)}
    resize!(ls, length(transitions))
    ls[1] = l1
    @inbounds for i in 2:length(transitions)
        t = transitions[i]
        ls[i] = convolve(build_line(t, params), instrument_function(t.wavenumber, params))
    end
    return ls
end
