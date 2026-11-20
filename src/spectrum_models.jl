using ComponentArrays

"""
    PhotonRate()

`amplitude` convention for [`EmissionSpectrum`](@ref) (default): bare
`A_ul · f_u`, a photon emission rate: the weighting of photon-counting
detectors (CCDs, PMTs).
"""
struct PhotonRate end

"""
    SpectralRadiance()

`amplitude` convention for [`EmissionSpectrum`](@ref): photon-energy weighted
`(h·c·ν̃/4π) · A_ul · f_u`, i.e. [`RadianceSpectrum`](@ref)'s emission
coefficient `j`. With `number_density_function`/`path_length_function` set, the
spectrum is the optically-thin limit of `RadianceSpectrum`, in
W·cm⁻²·sr⁻¹·(cm⁻¹)⁻¹.
"""
struct SpectralRadiance end

# `true` leaves the type of the amplitude unchanged.
@inline _amplitude_factor(::PhotonRate, ν̃) = true
@inline _amplitude_factor(::SpectralRadiance, ν̃) = h * c0_cm * ν̃ / 4π

add_background!(_, ::Nothing, args...) = nothing
add_background!(spectrum, bg::Number, args...) = spectrum .+= bg
function add_background!(spectrum, bg, xs, params)
    return spectrum .+= bg(xs, params)
end

_background_value(::Nothing, xs, params) = 0.0
_background_value(bg::Number, xs, params) = bg
_background_value(bg, xs, params) = bg(xs, params)

"""
    StaticLineshapes(f)

Marks a `params`-independent lineshape function: `f(t) -> LineShape` takes only
the transition, instead of the default `(t, params) -> LineShape`.

```julia
SpeciesModel(; transitions, distribution_function,
               lineshape_function = StaticLineshapes(t -> Voigt(σ_D(t), γ)), ...)
```

The models then build the shapes and their instrument convolutions once and
update only the amplitudes on later calls, which speeds up fitting loops, also
under `ForwardDiff`.

`transitions` must not be mutated in place. An `instrument_function` that
depends on `params` is allowed but loses the speed-up: whenever the first
line's convolved shape changes, all shapes are rebuilt.
"""
struct StaticLineshapes{F}
    f::F
end

@inline _shape(lineshape_function, t, params) = lineshape_function(t, params)
@inline _shape(lineshape_function::StaticLineshapes, t, params) = lineshape_function.f(t)

# Line centre: the transition wavenumber plus the species' `lineshift_function`.
@inline _position(t, ::Nothing, params) = t.wavenumber
@inline _position(t, lineshift_function, params) = t.wavenumber + lineshift_function(t, params)

@inline function _emission_amplitude(t, distribution_function, params::ComponentVector{T}, amplitude) where {T}
    f_u = max(zero(T), distribution_function(t.upper, params))
    return _amplitude_factor(amplitude, t.wavenumber) * t.A * f_u
end

@inline function _absorption_amplitude(t, distribution_function, params::ComponentVector{T}) where {T}
    # Not clamped at zero: absorbance depends only on the population difference
    # `f_l - f_u*g_l/g_u`, so a signed `f_l` carries gain (inverted lines).
    # Emission does clamp: a negative photon rate has no interpretation.
    f_l = distribution_function(t.lower, params)
    f_u = distribution_function(t.upper, params)
    return absorption_cross_section(t) * (f_l - f_u * degeneracy(t.lower) / degeneracy(t.upper))
end

"""
    emission_line(t, distribution_function, lineshape_function, params,
                            amplitude=PhotonRate())

Uninstrumented emission line of transition `t`: amplitude `A_ul · f_u`, with
`f_u = distribution_function(t.upper, params)` the upper-state population
(fraction or number density), times `h·c·ν̃/4π` when `amplitude` is
[`SpectralRadiance`](@ref). `f_u` is clamped at zero.
"""
@inline function emission_line(t, distribution_function, lineshape_function,
                                         params::ComponentVector{T}, amplitude=PhotonRate();
                                         position=t.wavenumber) where {T}
    amp = _emission_amplitude(t, distribution_function, params, amplitude)
    return Line(one(T) * position, one(T) * amp, _shape(lineshape_function, t, params))
end

"""
    absorption_line(t, distribution_function, lineshape_function, params)

Uninstrumented absorption line of transition `t`: amplitude
`σ₀ · (f_l − f_u·g_l/g_u)`, the integrated cross section
([`absorption_cross_section`](@ref)) weighted by the population difference, net
of stimulated emission.

`f_l` and `f_u` are not clamped at zero: only their difference is identifiable,
and a negative value expresses gain (inverted lines).
"""
@inline function absorption_line(t, distribution_function, lineshape_function,
                                           params::ComponentVector{T}; position=t.wavenumber) where {T}
    amp = _absorption_amplitude(t, distribution_function, params)
    return Line(one(T) * position, one(T) * amp, _shape(lineshape_function, t, params))
end

"""
    planck_radiance(ν̃, T) -> Float64

Planck's law: spectral radiance of a blackbody at `T` [K] versus wavenumber `ν̃`
[cm⁻¹], in W·cm⁻²·sr⁻¹·(cm⁻¹)⁻¹, the units of `j(ν)/α(ν)` in
[`RadianceSpectrum`](@ref).
"""
@inline planck_radiance(ν̃, T) = 2 * h * c0_cm^2 * ν̃^3 / expm1(h * c0_cm * ν̃ / (kB * T))

# Like `Line`, but with two amplitudes (α and j) sharing one position and
# shape, so RadianceSpectrum's non-LTE path evaluates the lineshape once per
# (transition, gridpoint) instead of twice.
struct _DualLine{T, LS<:LineShape{T}}
    position::T
    amp_α::T
    amp_j::T
    shape::LS
end

function _DualLine(position, amp_α, amp_j, shape::LineShape{T}) where {T}
    num_type = promote_type(typeof(position), typeof(amp_α), typeof(amp_j), T)
    new_shape = numtype_convert(num_type, shape)
    return _DualLine{num_type, typeof(new_shape)}(position, amp_α, amp_j, new_shape)
end

convolve(dl::_DualLine, ls::LineShape) = _DualLine(dl.position, dl.amp_α, dl.amp_j, convolve(dl.shape, ls))

"""
    radiance_line(t, distribution_function, lineshape_function, params)

Line of transition `t` carrying two amplitudes on one shape: the absorption
amplitude of [`absorption_line`](@ref) plus the emission amplitude
`(h·c·ν̃/4π)·A·f_u`. The photon-energy and per-steradian factors are always
applied here: `j` must match [`planck_radiance`](@ref)'s units for `j/α` to be a valid
source function.

Populations are clamped at zero, unlike [`absorption_line`](@ref).
"""
@inline function radiance_line(t, distribution_function, lineshape_function,
                                     params::ComponentVector{T}; position=t.wavenumber) where {T}
    f_l = max(zero(T), distribution_function(t.lower, params))
    f_u = max(zero(T), distribution_function(t.upper, params))
    g_l = degeneracy(t.lower)
    g_u = degeneracy(t.upper)
    amp_α = absorption_cross_section(t) * (f_l - f_u * g_l / g_u)
    amp_j = (h * c0_cm * t.wavenumber / 4π) * t.A * f_u
    return _DualLine(one(T) * position, one(T) * amp_α, one(T) * amp_j, _shape(lineshape_function, t, params))
end

@inline function _add_dual_line!(αret, jret, xs, dl::_DualLine, lo::Int, hi::Int)
    @inbounds for i in lo:hi
        v = dl.shape(xs[i] - dl.position)
        αret[i] += dl.amp_α * v
        jret[i] += dl.amp_j * v
    end
    return nothing
end

@inline function _add_dual_line_binned!(αret, jret, xs, dl::_DualLine, lo::Int, hi::Int, half_Δx, inv_Δx)
    # `BinAverage` returns a bin average, so both amplitudes carry 1/bin_width.
    ampα = dl.amp_α * inv_Δx
    ampj = dl.amp_j * inv_Δx
    @inbounds for i in lo:hi
        Δ = xs[i] - dl.position
        v = integrate(dl.shape, Δ - half_Δx, Δ + half_Δx)
        αret[i] += ampα * v
        jret[i] += ampj * v
    end
    return nothing
end

"""
    calculate_radiance_lines!(αret, jret, xs, dual_lines, method::Union{PointSample,BinAverage})

In-place [`calculate_radiance_lines`](@ref): fills `αret`/`jret` with the
absorption/emission coefficient sums over `dual_lines`.

Only `PointSample` and `BinAverage` are supported; there is no LDM equivalent for
dual-amplitude (α, j) lines. Their `min_amplitude` must be unset: it has no single
meaning for lines with two amplitudes. `αret` and `jret` must share indices with
`xs`, and `wing_cut` needs an ascending grid.
"""
function calculate_radiance_lines!(αret, jret, xs, dual_lines, method::PointSample)
    _reject_min_amplitude(method)
    _check_radiance_grid(αret, jret, xs, method)
    fill!(αret, zero(eltype(αret)))
    fill!(jret, zero(eltype(jret)))
    for dl in dual_lines
        lo, hi = line_index_range(xs, dl, method.wing_cut)
        lo > hi && continue
        _add_dual_line!(αret, jret, xs, dl, lo, hi)
    end
    return nothing
end

function calculate_radiance_lines!(αret, jret, xs, dual_lines, method::BinAverage)
    _reject_min_amplitude(method)
    _check_radiance_grid(αret, jret, xs, method)
    fill!(αret, zero(eltype(αret)))
    fill!(jret, zero(eltype(jret)))
    half_Δx = method.bin_width / 2
    inv_Δx = inv(method.bin_width)
    for dl in dual_lines
        lo, hi = line_index_range(xs, dl, method.wing_cut)
        lo > hi && continue
        _add_dual_line_binned!(αret, jret, xs, dl, lo, hi, half_Δx, inv_Δx)
    end
    return nothing
end

# The loops write both buffers under `@inbounds`, and `wing_cut` searches `xs`.
function _check_radiance_grid(αret, jret, xs, method)
    check_axes(αret, xs)
    check_axes(jret, xs)
    method.wing_cut === nothing || check_ascending(xs, "calculate_radiance_lines! with wing_cut")
    return nothing
end

_reject_min_amplitude(method) = method.min_amplitude === nothing || throw(ArgumentError(
    "`min_amplitude` is not supported on RadianceSpectrum's non-LTE path: each line " *
    "has an absorption and an emission amplitude in different units. Filter the " *
    "transitions instead."))

calculate_radiance_lines!(αret, jret, xs, dual_lines, method::LineSynthesisMethod) = error(
    "calculate_radiance_lines!(..., method=$(typeof(method))) is not supported: " *
    "only `PointSample` and `BinAverage` synthesize the dual-amplitude (α, j) lines of " *
    "RadianceSpectrum's non-LTE path; there is no LDM equivalent. " *
    "The LTE fast path does support any method."
)

"""
    calculate_radiance_lines(xs, dual_lines, method::Union{PointSample,BinAverage}) -> (α, j)

Allocating [`calculate_radiance_lines!`](@ref). Returns the absorption
coefficient `α` and emission coefficient `j`, before multiplication by number
density and path length (see [`RadianceSpectrum`](@ref)).
"""
function calculate_radiance_lines(xs, dual_lines, method::LineSynthesisMethod)
    isempty(dual_lines) && throw(ArgumentError(
        "calculate_radiance_lines: `dual_lines` is empty: need at least one line"))
    T = typeof(dual_lines[1].amp_α)
    αret = zeros(T, length(xs))
    jret = zeros(T, length(xs))
    calculate_radiance_lines!(αret, jret, xs, dual_lines, method)
    return αret, jret
end

# (1 − e^{−τ})/τ, by its series below |τ| = 1e-4. The series keeps the value
# and its ForwardDiff partials finite at τ = 0.
@inline function _weight_one_minus_exp(τ)
    abs(τ) < 1e-4 && return evalpoly(τ, (1.0, -0.5, inv(6.0), -inv(24.0)))
    return -expm1(-τ) / τ
end

# One grid point of the RTE:  I = (j/α)·(1 − e^{−αL}) + I₀·e^{−αL}, with the
# emission term written as j·L·(1 − e^{−τ})/τ so j/α is never formed (0/0).
@inline function _combine_radiance(j, α, L, I₀)
    τ = α * L
    return j * L * _weight_one_minus_exp(τ) + I₀ * exp(-τ)
end

# Models nonlinear in α (RadianceSpectrum, TransmittanceSpectrum) pass their
# instrument_function through this check.
_reject_instrumental(::Nothing) = nothing
_reject_instrumental(shape) = throw(ArgumentError(
    "instrument_function returned $(typeof(shape)). Per-line instrumental " *
    "broadening would convolve α before this model's nonlinear step in α, " *
    "which is wrong once the slab is not optically thin. Leave " *
    "instrument_function at its default and apply the instrument profile " *
    "to the model output with convolve_spectrum instead."))

# Per-species synthesis workspace, so FFT buffers and plans persist across
# calls. `nothing` for methods without one (PointSample, BinAverage).
synthesis_workspace(::LineSynthesisMethod) = nothing
synthesis_workspace(::LDM) = LDMWorkspace()
synthesis_workspace(::SparseLDM) = SparseLDMWorkspace()

_calculate_lines!(dest, ν, ls, method, ::Nothing) = calculate_lines!(dest, ν, ls, method)
_calculate_lines!(dest, ν, ls, method, ws) = calculate_lines!(dest, ν, ls, method, ws)

"""
    SpeciesModel(; transitions, distribution_function, lineshape_function,
                   number_density_function=(p) -> 1.0, lineshift_function=nothing,
                   method=PointSample())

One radiating/absorbing species: its line list, population distribution,
lineshape, number density [cm⁻³], and line-synthesis `method`. Every spectrum
model takes a `Tuple` of these via its `species` keyword; the coefficients add
before anything nonlinear, so several species share one slab:

```julia
sp = RadianceSpectrum(;
    species = (
        SpeciesModel(; transitions = tl_OH, distribution_function = bd_OH,
                       lineshape_function = lsf_OH,
                       number_density_function = p -> p.N_OH),
        SpeciesModel(; transitions = tl_H2O, distribution_function = bd_H2O,
                       lineshape_function = lsf_H2O,
                       number_density_function = p -> p.N_H2O,
                       method = LDM()),
    ),
    path_length_function = p -> p.L,
)
```

`method` is per-species, so a sparse atomic system (`PointSample`) and a dense
molecular band (`LDM`/`SparseLDM`) can share a slab. For `LDM`/`SparseLDM` a
synthesis workspace is created automatically and reused across calls. A model's
flat per-species keywords are shorthand for one `SpeciesModel` with the default
`PointSample()`.

`lineshift_function`: `(t, params) -> Δν̃` [cm⁻¹], added to every line
centre; `nothing` (default) leaves the centres at `t.wavenumber`. The pressure
shift of a HITRAN line list is `(t, p) -> t.pressure_shift * p.P` with `p.P`
in atm. A shift that depends on `params` disables the [`StaticLineshapes`](@ref)
fast path, since line positions then change between calls.

Slab and instrument properties (`path_length_function`, `incident_radiance`,
`background`, `wavenumber_function`, `instrument_function`, `amplitude`) stay
on the model.
"""
Base.@kwdef struct SpeciesModel{TR, DF, LF, NF, SF, M<:LineSynthesisMethod, W}
    transitions::Vector{TR}
    distribution_function::DF
    lineshape_function::LF
    number_density_function::NF = (p) -> 1.0
    lineshift_function::SF = nothing
    method::M = PointSample()
    # Reused across calls; not thread-safe.
    workspace::W = synthesis_workspace(method)
end

Base.show(io::IO, s::SpeciesModel) =
    print(io, "SpeciesModel(", length(s.transitions), length(s.transitions) == 1 ? " transition)" : " transitions)")

function _as_species_tuple(model::Symbol, species)
    species isa SpeciesModel && return (species,)
    sp = Tuple(species)
    isempty(sp) && throw(ArgumentError("$model: `species` is empty: need at least one SpeciesModel"))
    all(s -> s isa SpeciesModel, sp) || throw(ArgumentError(
        "$model: `species` must contain only SpeciesModel entries"))
    return sp
end

# Turn either a `species` tuple or the flat single-species keywords into a
# `Tuple{Vararg{SpeciesModel}}`.
function _normalize_species(model::Symbol, species, transitions, distribution_function,
                            lineshape_function, number_density_function)
    if species !== nothing
        if !(transitions === nothing && distribution_function === nothing &&
             lineshape_function === nothing && number_density_function === nothing)
            throw(ArgumentError("$model: pass either `species` or the flat per-species " *
                                "keywords (transitions/distribution_function/lineshape_function/" *
                                "number_density_function), not both"))
        end
        return _as_species_tuple(model, species)
    end
    transitions === nothing && throw(ArgumentError("$model: `transitions` (or `species`) is required"))
    distribution_function === nothing && throw(ArgumentError("$model: `distribution_function` is required"))
    lineshape_function === nothing && throw(ArgumentError("$model: `lineshape_function` is required"))
    return (SpeciesModel(; transitions, distribution_function, lineshape_function,
                         number_density_function = something(number_density_function, (p) -> 1.0)),)
end

"""
    Slab(; species, path_length_function)

One uniform slab of a line of sight: a `Tuple` of [`SpeciesModel`](@ref) (or a
single one) plus its geometric length `(params) -> L` [cm]. The geometric
models take a `Tuple` of these via their `slabs` keyword; `slabs[1]` is the
slab farthest from the observer, and light exits `slabs[end]`:

```julia
sp = RadianceSpectrum(; slabs = (
    Slab(; species = core_species, path_length_function = p -> p.L_core),  # back
    Slab(; species = edge_species, path_length_function = p -> p.L_edge),  # front
))
```

Per-slab conditions come from the parameter keys the species' functions capture
(`p.T_core` vs `p.T_edge` via `Boltzmann`'s `temperature_function`, `p.N_core`
vs `p.N_edge`, …). The `species` + `path_length_function` model keywords are
shorthand for a one-slab tuple.
"""
struct Slab{SP<:Tuple, PF}
    species::SP
    path_length_function::PF
end

Slab(; species, path_length_function) =
    Slab(_as_species_tuple(:Slab, species), path_length_function)

Base.show(io::IO, s::Slab) =
    print(io, "Slab(", length(s.species), " species)")

# Normalize to `Tuple{Vararg{Slab}}`: either `slabs` directly, or the
# `species`/flat keywords plus `path_length_function` as a one-slab shorthand.
function _normalize_slabs(model::Symbol, slabs, species, transitions, distribution_function,
                          lineshape_function, number_density_function,
                          path_length_function, default_L)
    if slabs !== nothing
        if !(species === nothing && transitions === nothing && distribution_function === nothing &&
             lineshape_function === nothing && number_density_function === nothing &&
             path_length_function === nothing)
            throw(ArgumentError("$model: pass either `slabs` or the per-slab keywords " *
                                "(`species`/flat keywords and `path_length_function`), not both"))
        end
        slabs isa Slab && return (slabs,)
        sl = Tuple(slabs)
        isempty(sl) && throw(ArgumentError("$model: `slabs` is empty: need at least one Slab"))
        all(s -> s isa Slab, sl) || throw(ArgumentError(
            "$model: `slabs` must contain only Slab entries"))
        return sl
    end
    sp = _normalize_species(model, species, transitions, distribution_function,
                            lineshape_function, number_density_function)
    L = path_length_function === nothing ? default_L : path_length_function
    L === nothing && throw(ArgumentError("$model: `path_length_function` is required (or pass `slabs`)"))
    return (Slab(sp, L),)
end

# The common temperature of an all-LTE species tuple, or `nothing`. The
# temperatures are compared by value at these params.
function _lte_temperature(species::Tuple, params)
    _all_lte(species) || return nothing
    T1 = _temperature(species[1].distribution_function, params)
    return _same_temperature(T1, Base.tail(species), params) ? T1 : nothing
end
_same_temperature(T1, ::Tuple{}, params) = true
_same_temperature(T1, species::Tuple, params) =
    _temperature(species[1].distribution_function, params) == T1 &&
    _same_temperature(T1, Base.tail(species), params)

# One `LineCache` per species. The caches belong to the model: a static build
# folds the model's instrument into the cached shapes.
_species_caches(species::Tuple) = map(_ -> LineCache(), species)
_slab_caches(slabs::Tuple) = map(slab -> _species_caches(slab.species), slabs)

function _fill_species_emission!(cache::LineCache, s::SpeciesModel, instrument_function, amplitude, params)
    # A line shift moves the centres between calls, so it takes the full path.
    # In the static build only the amplitude carries the params eltype.
    if s.lineshape_function isa StaticLineshapes && s.lineshift_function === nothing
        lsf = s.lineshape_function
        return update_or_fill_lines!(cache, :emission, :emission_amps, s.transitions,
                    (t, p) -> Line(t.wavenumber, _emission_amplitude(t, s.distribution_function, p, amplitude), lsf.f(t)),
                    (t, p) -> _emission_amplitude(t, s.distribution_function, p, amplitude),
                    instrument_function, params)
    end
    return fill_lines!(cache, :emission, s.transitions,
                (t, p) -> emission_line(t, s.distribution_function, s.lineshape_function, p, amplitude;
                                                  position = _position(t, s.lineshift_function, p)),
                instrument_function, params)
end

function _fill_species_absorption!(cache::LineCache, s::SpeciesModel, instrument_function, params)
    if s.lineshape_function isa StaticLineshapes && s.lineshift_function === nothing
        lsf = s.lineshape_function
        return update_or_fill_lines!(cache, :absorption, :absorption_amps, s.transitions,
                    (t, p) -> Line(t.wavenumber, _absorption_amplitude(t, s.distribution_function, p), lsf.f(t)),
                    (t, p) -> _absorption_amplitude(t, s.distribution_function, p),
                    instrument_function, params)
    end
    return fill_lines!(cache, :absorption, s.transitions,
                (t, p) -> absorption_line(t, s.distribution_function, s.lineshape_function, p;
                                                    position = _position(t, s.lineshift_function, p)),
                instrument_function, params)
end

_fill_species_dual!(cache::LineCache, s::SpeciesModel, instrument_function, params) =
    fill_lines!(cache, :dual, s.transitions,
                (t, p) -> radiance_line(t, s.distribution_function, s.lineshape_function, p;
                                              position = _position(t, s.lineshift_function, p)),
                instrument_function, params)

_all_lte(species::Tuple) = all(s -> is_LTE(s.distribution_function), species)

# Output eltype for the allocating forms: the promotion of the grid, every
# species' lines and every species' number density.
_species_sum_eltype(ν, species, line_vectors, params) =
    promote_type(map(ls -> spectrum_eltype(ν, ls), line_vectors)...,
                 map(s -> typeof(s.number_density_function(params)), species)...)

# Σₛ Nₛ·(line sum of species s) into `dest`.
function _sum_species!(dest, cache::LineCache, ν, species::Tuple, line_vectors::Tuple, params)
    if length(species) === 1
        _calculate_lines!(dest, ν, line_vectors[1], species[1].method, species[1].workspace)
        N = species[1].number_density_function(params)
        isone(N) || (dest .*= N)
    else
        fill!(dest, zero(eltype(dest)))
        tmp = _scratch_vector!(cache, :species_tmp, zero(eltype(dest)), length(ν))
        foreach(species, line_vectors) do s, ls
            _calculate_lines!(tmp, ν, ls, s.method, s.workspace)
            N = s.number_density_function(params)
            @. dest += N * tmp
        end
    end
    return dest
end

# Output eltype for the multi-slab allocating forms, found without building
# lines: the promotion of the grid, params, N and L types.
_slabs_eltype(ν, slabs::Tuple, params) =
    promote_type(Float64, eltype(params), typeof(first(ν)),
                 map(slab -> promote_type(typeof(slab.path_length_function(params)),
                                          map(s -> typeof(s.number_density_function(params)), slab.species)...),
                     slabs)...)

# Σₖ Lₖ · (species sum of slab k) into `dest`, for the models that sum lines
# over slabs.
function _sum_slabs!(dest, cache::LineCache, ν, slabs::Tuple, caches::Tuple, fill_species, params)
    if length(slabs) === 1
        slab = slabs[1]
        line_vectors = map(fill_species, slab.species, caches[1])
        _sum_species!(dest, cache, ν, slab.species, line_vectors, params)
        L = slab.path_length_function(params)
        isone(L) || (dest .*= L)
    else
        fill!(dest, zero(eltype(dest)))
        tmp = _scratch_vector!(cache, :slab_tmp, zero(eltype(dest)), length(ν))
        foreach(slabs, caches) do slab, species_caches
            line_vectors = map(fill_species, slab.species, species_caches)
            _sum_species!(tmp, cache, ν, slab.species, line_vectors, params)
            L = slab.path_length_function(params)
            @. dest += L * tmp
        end
    end
    return dest
end

# The dual (α, j) analogue for the non-LTE radiance path.
function _sum_dual_species!(αdest, jdest, cache::LineCache, ν, species::Tuple, dual_line_vectors::Tuple, params)
    if length(species) === 1
        calculate_radiance_lines!(αdest, jdest, ν, dual_line_vectors[1], species[1].method)
        N = species[1].number_density_function(params)
        if !isone(N)
            αdest .*= N
            jdest .*= N
        end
    else
        fill!(αdest, zero(eltype(αdest)))
        fill!(jdest, zero(eltype(jdest)))
        αtmp = _scratch_vector!(cache, :α_species_tmp, zero(eltype(αdest)), length(ν))
        jtmp = _scratch_vector!(cache, :j_species_tmp, zero(eltype(jdest)), length(ν))
        foreach(species, dual_line_vectors) do s, dls
            calculate_radiance_lines!(αtmp, jtmp, ν, dls, s.method)
            N = s.number_density_function(params)
            @. αdest += N * αtmp
            @. jdest += N * jtmp
        end
    end
    return αdest, jdest
end

"""
    EmissionSpectrum(; species, path_length_function=(p) -> 1.0,
                       instrument_function=(wn, p) -> nothing,
                       wavenumber_function=(x, p) -> x,
                       amplitude=PhotonRate(),
                       background=nothing)
    EmissionSpectrum(; slabs, ...)

Optically-thin emission model: `Σₖ Lₖ · Σₛ Nₛ · Σ A_ul · f_u(state) · φ(ν)` over
the transitions of every [`SpeciesModel`](@ref), with no self-absorption (so
slabs add). Call as `sp(xs, params)` (allocating) or `sp(ret, xs, params)`
(in-place), with `params` a `ComponentVector`. [`RadianceSpectrum`](@ref) is the
self-absorbing counterpart.

Instead of `species`, the flat keywords `transitions`, `distribution_function`,
`lineshape_function`, `number_density_function` build a one-species tuple with
the default `PointSample()` method; for another method, pass a `SpeciesModel`.

# Keyword arguments

- `species`: `Tuple` of [`SpeciesModel`](@ref) (or one `SpeciesModel`)
- `path_length_function`: `(params) -> L` [cm], default `1.0`. Shared by all
  species; each species carries its own number density.
- `slabs`: `Tuple` of [`Slab`](@ref) for a multi-slab line of sight, mutually
  exclusive with `species`/`path_length_function` (one-slab sugar).
- `instrument_function`: `(wavenumber, params) -> LineShape | nothing`;
  `nothing` (default) applies no instrumental broadening. Instrument profiles
  with no closed-form combination with the lineshape (e.g. `Triangular`) go
  through [`convolve_spectrum`](@ref) instead.
- `wavenumber_function`: `(xs, params) -> wavenumbers`, for wavelength
  calibration or axis remapping. Default identity.
- `amplitude`: [`PhotonRate()`](@ref PhotonRate) (default) or
  [`SpectralRadiance()`](@ref SpectralRadiance). Across a wide band the two
  differ in relative line intensities (a factor `ν̃`), not just units.
- `background`: constant or `(xs, params) -> offset`, added to the spectrum

# Example

```julia
# CO (HITRAN molecule 5, isotopologue 1) from a local .par file
lines = HITRAN.read_par("05_HITEMP2019.par"; wn_min = 2100, wn_max = 2200)
tl = make_transition_list(filter(l -> l.iso == 1, lines))
m = molecular_mass_kg(HITRAN.molparam(5, 1))
sp = EmissionSpectrum(;
    transitions = tl,
    distribution_function = Boltzmann(HITRAN.partition_function(5, 1), 0.0, State(:X)),
    lineshape_function = (t, p) -> Voigt(sigma_doppler(t.wavenumber, p.T; m), p.gamma),
    instrument_function = (wn, p) -> Gaussian(p.sigma_inst),
)
spectrum = sp(2100:0.01:2200, ComponentVector(T = 1500.0, gamma = 0.05, sigma_inst = 0.1))
```
"""
struct EmissionSpectrum{SL<:Tuple, IF, WF, AM, B, SC<:Tuple}
    slabs::SL
    instrument_function::IF
    wavenumber_function::WF
    amplitude::AM
    line_cache::LineCache   # grid scratch slots
    background::B
    species_caches::SC      # per-species line vectors
end

function EmissionSpectrum(; slabs = nothing, species = nothing,
                            transitions = nothing, distribution_function = nothing,
                            lineshape_function = nothing, number_density_function = nothing,
                            instrument_function = (wn, p) -> nothing,
                            wavenumber_function = (x, p) -> x,
                            path_length_function = nothing,
                            amplitude = PhotonRate(),
                            line_cache = LineCache(),
                            background = nothing)
    sl = _normalize_slabs(:EmissionSpectrum, slabs, species, transitions, distribution_function,
                          lineshape_function, number_density_function,
                          path_length_function, (p) -> 1.0)
    return EmissionSpectrum(sl, instrument_function, wavenumber_function,
                            amplitude, line_cache, background, _slab_caches(sl))
end

function (sp::EmissionSpectrum)(xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    spectrum = zeros(_slabs_eltype(ν, sp.slabs, params), length(ν))
    _sum_slabs!(spectrum, sp.line_cache, ν, sp.slabs, sp.species_caches,
                           (s, c) -> _fill_species_emission!(c, s, sp.instrument_function, sp.amplitude, params), params)
    add_background!(spectrum, sp.background, xs, params)
    return spectrum
end

function (sp::EmissionSpectrum)(ret, xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    _sum_slabs!(ret, sp.line_cache, ν, sp.slabs, sp.species_caches,
                           (s, c) -> _fill_species_emission!(c, s, sp.instrument_function, sp.amplitude, params), params)
    add_background!(ret, sp.background, xs, params)
    return ret
end

"""
    RadianceSpectrum(; species, path_length_function,
                        instrument_function=(wn, p) -> nothing,
                        wavenumber_function=(x, p) -> x,
                        incident_radiance=nothing, background=nothing)
    RadianceSpectrum(; slabs, ...)

Self-absorbing radiative transfer through uniform, non-scattering slabs. Per
slab, computes `α(ν) = Σₛ αₛ(ν)` and `j(ν) = Σₛ jₛ(ν)` over its
[`SpeciesModel`](@ref)s and folds
`I(ν) ← S(ν)·(1 − e^{−α(ν)L}) + I(ν)·e^{−α(ν)L}` with source function
`S(ν) = j(ν)/α(ν)`, starting from `incident_radiance` behind `slabs[1]` and
exiting `slabs[end]` toward the observer. Coefficients add before the RTE, so
overlapping lines of several species self-absorb against each other.
[`EmissionSpectrum`](@ref) is the optically-thin
counterpart.

# Keyword arguments

`species` (or the flat single-species keywords, see [`SpeciesModel`](@ref)) and
`wavenumber_function` as for [`EmissionSpectrum`](@ref), plus:

- `slabs`: `Tuple` of [`Slab`](@ref), ordered far side first; mutually
  exclusive with `species`/`path_length_function` (one-slab sugar).
- `path_length_function`: `(params) -> L` [cm]; required in the one-slab form.
  Each species carries its own `number_density_function`, e.g.
  `(p) -> p.P * p.x / (kB * p.T) / 1e6` for pressure `p.P` [Pa], mole fraction
  `p.x`, temperature `p.T` [K].
- `instrument_function`: must return `nothing` (the default); a returned
  shape throws. Use [`convolve_spectrum`](@ref), which
  convolves after the radiative transfer.
- `incident_radiance`: `I₀(ν)` entering behind `slabs[1]`, attenuated by
  `e^{−αL}`; constant, callable `(ν, params) -> ...` on the mapped wavenumber
  axis, or `nothing` (default, `I₀ = 0`)
- `background`: instrumental offset (dark signal, straylight) added to the
  final spectrum on the raw `xs` axis. Not attenuated by the slab.

Per-species `method` is unrestricted only on the LTE fast path below; the
non-LTE path supports `PointSample`/`BinAverage` (see
[`calculate_radiance_lines!`](@ref)).

With `BinAverage`, `α` and `j` are bin-averaged before the radiative transfer. That is
not the bin-averaged radiance; the error is second order in the bin width.

# LTE fast path

Decided per slab: when every species' distribution satisfies `is_LTE` (an
unwrapped `Boltzmann`) and their temperatures are equal at the current
`params`, the source function is [`planck_radiance`](@ref) in closed form and `j` is
never computed. Otherwise
the slab takes the general (α, j) path. The
two paths differ at the 1e-4 level across a line: Planck is evaluated at each
grid point here, and at the line centre in the general path.

A bare `Boltzmann` whose `parent` excludes lower states of its species'
transitions throws at construction: those populations are zero, which makes the
absorption coefficient negative.
"""
struct RadianceSpectrum{SL<:Tuple, IF, WF, IR, B, SC<:Tuple}
    slabs::SL
    instrument_function::IF
    wavenumber_function::WF
    incident_radiance::IR
    background::B
    line_cache::LineCache   # grid scratch slots
    species_caches::SC      # per-species line vectors
end

function RadianceSpectrum(; slabs = nothing, species = nothing,
                            transitions = nothing, distribution_function = nothing,
                            lineshape_function = nothing, number_density_function = nothing,
                            path_length_function = nothing,
                            instrument_function = (wn, p) -> nothing,
                            wavenumber_function = (x, p) -> x,
                            incident_radiance = nothing,
                            background = nothing,
                            line_cache = LineCache())
    sl = _normalize_slabs(:RadianceSpectrum, slabs, species, transitions, distribution_function,
                          lineshape_function, number_density_function,
                          path_length_function, nothing)
    foreach(slab -> foreach(_check_boltzmann_parent, slab.species), sl)
    return RadianceSpectrum(sl, instrument_function, wavenumber_function,
                            incident_radiance, background, line_cache, _slab_caches(sl))
end

# A bare `Boltzmann` gives zero population outside its `parent`. A lower state
# out there makes α negative, so it is rejected. An upper state outside only
# drops stimulated emission.
function _check_boltzmann_parent(s::SpeciesModel)
    df = s.distribution_function
    df isa Boltzmann || return nothing
    all(t -> t.lower in df.parent, s.transitions) || throw(ArgumentError(
        "RadianceSpectrum: the Boltzmann parent $(df.parent) excludes lower states of " *
        "this species' transitions, so their population is zero and the absorption " *
        "coefficient negative. Use a distribution that covers the lower states too, " *
        "or EmissionSpectrum for the optically thin case."))
    return nothing
end

_reject_wrapped(sp) = (wn, p) -> _reject_instrumental(sp.instrument_function(wn, p))

# One slab step of the RTE fold: I ← S·(1 − e^{−αL}) + I·e^{−αL}, elementwise.
function _apply_slab!(I, cache::LineCache, ν, slab::Slab, species_caches::Tuple, instr, params)
    L = slab.path_length_function(params)
    T_lte = _lte_temperature(slab.species, params)
    if T_lte !== nothing
        line_vectors = map((s, c) -> _fill_species_absorption!(c, s, instr, params), slab.species, species_caches)
        α = _scratch_vector!(cache, :α_scratch,
                             zero(_species_sum_eltype(ν, slab.species, line_vectors, params)), length(ν))
        _sum_species!(α, cache, ν, slab.species, line_vectors, params)
        @. I = planck_radiance(ν, T_lte) * -expm1(-α * L) + I * exp(-α * L)
    else
        dual_line_vectors = map((s, c) -> _fill_species_dual!(c, s, instr, params), slab.species, species_caches)
        S = promote_type(map(dls -> promote_type(typeof(first(ν)), typeof(dls[1].position), typeof(dls[1].amp_α)), dual_line_vectors)...,
                         map(s -> typeof(s.number_density_function(params)), slab.species)...)
        α = _scratch_vector!(cache, :α_scratch, zero(S), length(ν))
        j = _scratch_vector!(cache, :j_scratch, zero(S), length(ν))
        _sum_dual_species!(α, j, cache, ν, slab.species, dual_line_vectors, params)
        @. I = _combine_radiance(j, α, L, I)
    end
    return I
end

_fold_slabs!(I, cache, ν, ::Tuple{}, ::Tuple{}, instr, params) = I
function _fold_slabs!(I, cache, ν, slabs::Tuple, caches::Tuple, instr, params)
    _apply_slab!(I, cache, ν, slabs[1], caches[1], instr, params)
    return _fold_slabs!(I, cache, ν, Base.tail(slabs), Base.tail(caches), instr, params)
end

function (sp::RadianceSpectrum)(ret, xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    # Incident radiance is physical and lives on the mapped wavenumber axis; the
    # instrumental background below is added on the raw xs axis.
    I₀ = _background_value(sp.incident_radiance, ν, params)
    ret .= zero(eltype(ret))
    ret .+= I₀
    _fold_slabs!(ret, sp.line_cache, ν, sp.slabs, sp.species_caches, _reject_wrapped(sp), params)
    add_background!(ret, sp.background, xs, params)
    return ret
end

function (sp::RadianceSpectrum)(xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    I = zeros(_slabs_eltype(ν, sp.slabs, params), length(ν))
    I .+= _background_value(sp.incident_radiance, ν, params)
    _fold_slabs!(I, sp.line_cache, ν, sp.slabs, sp.species_caches, _reject_wrapped(sp), params)
    add_background!(I, sp.background, xs, params)
    return I
end

"""
    AbsorptionCoefficientSpectrum(; species,
                                    instrument_function=(wn, p) -> nothing,
                                    wavenumber_function=(x, p) -> x,
                                    background=nothing)

The absorption coefficient `α(ν) = Σₛ Nₛ · Σ σ₀·(f_l − f_u·g_l/g_u)·φ(ν)`
[cm⁻¹] over the [`SpeciesModel`](@ref)s in `species` (see
[`absorption_line`](@ref)). The building block behind
[`AbsorbanceSpectrum`](@ref), [`TransmittanceSpectrum`](@ref) and
[`RadianceSpectrum`](@ref)'s LTE path. The flat single-species keywords work as
for [`EmissionSpectrum`](@ref).

With a single species whose `number_density_function` is left at `1.0`, the
result is the effective cross-section spectrum [cm²] per molecule.

`α` is linear in the lineshape, so unlike the nonlinear models:

- `instrument_function` is allowed: per-line convolution is exact here;
- any per-species `method` works, including `LDM`;
- `BinAverage` returns the exact bin average.

`background` models continuum absorption, added to `α` on the raw `xs` axis.
"""
struct AbsorptionCoefficientSpectrum{SP<:Tuple, IF, WF, B, SC<:Tuple}
    species::SP
    instrument_function::IF
    wavenumber_function::WF
    background::B
    line_cache::LineCache   # grid scratch slots
    species_caches::SC      # per-species line vectors
end

function AbsorptionCoefficientSpectrum(; species = nothing,
                                         transitions = nothing, distribution_function = nothing,
                                         lineshape_function = nothing, number_density_function = nothing,
                                         instrument_function = (wn, p) -> nothing,
                                         wavenumber_function = (x, p) -> x,
                                         background = nothing,
                                         line_cache = LineCache())
    sp = _normalize_species(:AbsorptionCoefficientSpectrum, species, transitions, distribution_function,
                            lineshape_function, number_density_function)
    return AbsorptionCoefficientSpectrum(sp, instrument_function, wavenumber_function,
                                         background, line_cache, _species_caches(sp))
end

function (sp::AbsorptionCoefficientSpectrum)(xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    line_vectors = map((s, c) -> _fill_species_absorption!(c, s, sp.instrument_function, params),
              sp.species, sp.species_caches)
    α = zeros(_species_sum_eltype(ν, sp.species, line_vectors, params), length(ν))
    _sum_species!(α, sp.line_cache, ν, sp.species, line_vectors, params)
    add_background!(α, sp.background, xs, params)
    return α
end

function (sp::AbsorptionCoefficientSpectrum)(ret, xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    line_vectors = map((s, c) -> _fill_species_absorption!(c, s, sp.instrument_function, params),
              sp.species, sp.species_caches)
    _sum_species!(ret, sp.line_cache, ν, sp.species, line_vectors, params)
    add_background!(ret, sp.background, xs, params)
    return ret
end

"""
    AbsorbanceSpectrum(; species, path_length_function,
                         instrument_function=(wn, p) -> nothing,
                         wavenumber_function=(x, p) -> x,
                         background=nothing)

Napierian absorbance `A(ν) = Σₖ αₖ(ν)·Lₖ = −ln(I/I₀)`, i.e. the optical depth of
the line of sight; decadic absorbance is `A / ln(10)`. `α(ν)` sums over each
slab's [`SpeciesModel`](@ref)s; pass `species` + `path_length_function` (one
slab) or a `slabs` `Tuple` of [`Slab`](@ref) (absorbances add). The flat
single-species keywords work as for [`EmissionSpectrum`](@ref).

Linear in `α`, so `instrument_function` and any per-species `method` are
allowed, as in [`AbsorptionCoefficientSpectrum`](@ref). A spectrometer measures
transmittance, though: instrument-convolved absorbance matches an
instrument-convolved measurement only where the slab is optically thin or the
lines are resolved. Otherwise fit [`TransmittanceSpectrum`](@ref) +
[`convolve_spectrum`](@ref).

`background` is an additive baseline on the raw `xs` axis.
"""
struct AbsorbanceSpectrum{SL<:Tuple, IF, WF, B, SC<:Tuple}
    slabs::SL
    instrument_function::IF
    wavenumber_function::WF
    background::B
    line_cache::LineCache   # grid scratch slots
    species_caches::SC      # per-species line vectors
end

function AbsorbanceSpectrum(; slabs = nothing, species = nothing,
                              transitions = nothing, distribution_function = nothing,
                              lineshape_function = nothing, number_density_function = nothing,
                              path_length_function = nothing,
                              instrument_function = (wn, p) -> nothing,
                              wavenumber_function = (x, p) -> x,
                              background = nothing,
                              line_cache = LineCache())
    sl = _normalize_slabs(:AbsorbanceSpectrum, slabs, species, transitions, distribution_function,
                          lineshape_function, number_density_function,
                          path_length_function, nothing)
    return AbsorbanceSpectrum(sl, instrument_function, wavenumber_function,
                              background, line_cache, _slab_caches(sl))
end

function (sp::AbsorbanceSpectrum)(xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    A = zeros(_slabs_eltype(ν, sp.slabs, params), length(ν))
    _sum_slabs!(A, sp.line_cache, ν, sp.slabs, sp.species_caches,
                           (s, c) -> _fill_species_absorption!(c, s, sp.instrument_function, params), params)
    add_background!(A, sp.background, xs, params)
    return A
end

function (sp::AbsorbanceSpectrum)(ret, xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    _sum_slabs!(ret, sp.line_cache, ν, sp.slabs, sp.species_caches,
                           (s, c) -> _fill_species_absorption!(c, s, sp.instrument_function, params), params)
    add_background!(ret, sp.background, xs, params)
    return ret
end

"""
    TransmittanceSpectrum(; species, path_length_function,
                            wavenumber_function=(x, p) -> x,
                            background=nothing)

Transmittance `T(ν) = e^{−Σₖ αₖ(ν)·Lₖ}`, the normalized `I/I₀` of a non-emitting
line of sight. `α(ν)` sums over each slab's [`SpeciesModel`](@ref)s, so
overlapping absorbers and consecutive slabs multiply their transmittances
exactly. Pass `species` + `path_length_function` (one slab) or a `slabs` `Tuple`
of [`Slab`](@ref). The flat single-species keywords work as for
[`EmissionSpectrum`](@ref).

Nonlinear in `α`, hence the same two caveats as `RadianceSpectrum`:
`instrument_function` must return `nothing` (use [`convolve_spectrum`](@ref)), and `BinAverage` combines
bin-averaged `α` rather than bin-averaging `T` (see `RadianceSpectrum`).

`background` is an additive instrumental offset (e.g. straylight) on the raw
`xs` axis.
"""
struct TransmittanceSpectrum{SL<:Tuple, IF, WF, B, SC<:Tuple}
    slabs::SL
    instrument_function::IF
    wavenumber_function::WF
    background::B
    line_cache::LineCache   # grid scratch slots
    species_caches::SC      # per-species line vectors
end

function TransmittanceSpectrum(; slabs = nothing, species = nothing,
                                 transitions = nothing, distribution_function = nothing,
                                 lineshape_function = nothing, number_density_function = nothing,
                                 path_length_function = nothing,
                                 instrument_function = (wn, p) -> nothing,
                                 wavenumber_function = (x, p) -> x,
                                 background = nothing,
                                 line_cache = LineCache())
    sl = _normalize_slabs(:TransmittanceSpectrum, slabs, species, transitions, distribution_function,
                          lineshape_function, number_density_function,
                          path_length_function, nothing)
    return TransmittanceSpectrum(sl, instrument_function, wavenumber_function,
                                 background, line_cache, _slab_caches(sl))
end

function (sp::TransmittanceSpectrum)(xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    τ = zeros(_slabs_eltype(ν, sp.slabs, params), length(ν))
    _sum_slabs!(τ, sp.line_cache, ν, sp.slabs, sp.species_caches,
                           (s, c) -> _fill_species_absorption!(c, s, _reject_wrapped(sp), params), params)
    @. τ = exp(-τ)
    add_background!(τ, sp.background, xs, params)
    return τ
end

function (sp::TransmittanceSpectrum)(ret, xs, params::ComponentVector)
    ν = sp.wavenumber_function(xs, params)
    _sum_slabs!(ret, sp.line_cache, ν, sp.slabs, sp.species_caches,
                           (s, c) -> _fill_species_absorption!(c, s, _reject_wrapped(sp), params), params)
    @. ret = exp(-ret)
    add_background!(ret, sp.background, xs, params)
    return ret
end

function Base.show(io::IO, sp::AbsorptionCoefficientSpectrum)
    n_sp = length(sp.species)
    n_tr = sum(s -> length(s.transitions), sp.species)
    print(io, "AbsorptionCoefficientSpectrum(", n_sp, " species, ",
          n_tr, n_tr == 1 ? " transition)" : " transitions)")
end

for T in (:EmissionSpectrum, :RadianceSpectrum, :AbsorbanceSpectrum, :TransmittanceSpectrum)
    @eval function Base.show(io::IO, sp::$T)
        n_sl = length(sp.slabs)
        n_sp = sum(sl -> length(sl.species), sp.slabs)
        n_tr = sum(sl -> sum(s -> length(s.transitions), sl.species), sp.slabs)
        print(io, $(string(T)), "(")
        n_sl == 1 || print(io, n_sl, " slabs, ")
        print(io, n_sp, " species, ", n_tr, n_tr == 1 ? " transition)" : " transitions)")
    end
end
