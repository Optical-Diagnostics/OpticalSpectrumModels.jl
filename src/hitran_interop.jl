"""
    make_transition_list(lines::Vector{HITRAN.Line}) -> Vector{ElectricDipoleTransition}

Convert parsed HITRAN/HITEMP lines into `ElectricDipoleTransition`s. Uses only
the direct HITRAN160 fields `wav`, `A`, `El`, `gp`, `gpp`, so LTE modelling works
without decoding the `globu`/`globl`/`locu`/`locl` quantum-number strings.

Electronic/vibrational/rotational fields are left as `:X`/`nothing`
placeholders; decoding those strings is a prerequisite for driving non-LTE
distributions such as `StateResolved` from HITRAN.

The Einstein `A` is abundance-free, unlike HITRAN's intensity `S`, so
`number_density_function` must return the number density of the isotopologue
itself. For natural abundance, scale by
`HITRAN.molparam(mol_id, iso).abundance`. Lines of several isotopologues give a
warning: filter to one per `SpeciesModel`.

Lines without an Einstein `A`, statistical weights or a lower-state energy
(unassigned lines) throw an `ArgumentError`, since they would turn the spectrum
into `NaN`. Drop them first with `filter(HITRAN.is_usable, lines)`; they are not
dropped here, so that `lines[i]` keeps matching the `i`-th transition.

The air pressure shift `δ_air` becomes the transition's `pressure_shift`
(0 where the line list leaves it blank). Apply it with a `lineshift_function`
on the `SpeciesModel`, e.g. `(t, p) -> t.pressure_shift * p.P` for `P` in atm.
"""
function make_transition_list(lines::Vector{HITRAN.Line})
    n_bad = count(!HITRAN.is_usable, lines)
    n_bad == 0 || throw(ArgumentError(
        "make_transition_list: $n_bad of $(length(lines)) lines lack an Einstein A, " *
        "statistical weights or a lower-state energy (unassigned lines). Drop them " *
        "first: `filter(HITRAN.is_usable, lines)`."))
    # One SpeciesModel applies one population and number density to all its
    # lines, which is wrong across isotopologues.
    isos = unique((l.mol_id, l.iso) for l in lines)
    length(isos) > 1 && @warn "make_transition_list: the lines span " *
        "$(length(isos)) isotopologues $(isos). A SpeciesModel applies one " *
        "distribution and number density to all of them; give each isotopologue " *
        "its own, e.g. `filter(l -> l.iso == 1, lines)`."
    map(lines) do l
        lower = State(:X, nothing, nothing, l.El, l.gpp)
        upper = State(:X, nothing, nothing, l.El + l.wav, l.gp)
        ElectricDipoleTransition(lower, upper, l.wav, l.A, isnan(l.Pshft) ? 0.0 : l.Pshft)
    end
end
