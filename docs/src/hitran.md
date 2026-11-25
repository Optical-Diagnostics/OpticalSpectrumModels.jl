# HITRAN

`OpticalSpectrumModels.HITRAN` is a self-contained submodule for reading HITRAN/HITEMP
line-list files (the fixed-width "HITRAN160" `.par` format), the bundled
isotopologue metadata table, and the TIPS-2025 partition sums. It does not
depend on the rest of the package; `make_transition_list` converts its parsed
lines into [`ElectricDipoleTransition`](@ref OpticalSpectrumModels.ElectricDipoleTransition)s.

A HITRAN line list drives the spectrum models through `make_transition_list`
and a `Boltzmann` distribution with the TIPS partition sum:

```julia
lines = HITRAN.read_par("CO2.par"; wn_min = 2300.0, wn_max = 2400.0)
lines = filter(l -> l.iso == 1, lines)    # one isotopologue per species (626)
tl = make_transition_list(lines)
bd = Boltzmann(HITRAN.partition_function(2, 1), 0.0, State(:X))
```

A `.par` file holds every isotopologue of the molecule; give each its own
`SpeciesModel` with its own partition sum. The Einstein `A` in the line list is
abundance-free, so the number density is that of the isotopologue itself.

The pressure shift of each line is carried as `pressure_shift` and applied with
`SpeciesModel(; ..., lineshift_function = (t, p) -> t.pressure_shift * p.P)`,
with `p.P` in atm.

```@docs
OpticalSpectrumModels.HITRAN
OpticalSpectrumModels.HITRAN.Line
OpticalSpectrumModels.HITRAN.read_par
OpticalSpectrumModels.HITRAN.is_usable
OpticalSpectrumModels.HITRAN.read_hitemp
OpticalSpectrumModels.HITRAN.molecule_name
OpticalSpectrumModels.HITRAN.molparam
OpticalSpectrumModels.HITRAN.MOLPARAM
OpticalSpectrumModels.HITRAN.partition_function
OpticalSpectrumModels.HITRAN.PartitionFunction
OpticalSpectrumModels.HITRAN.TIPS
OpticalSpectrumModels.HITRAN.linestrength
```
