```@meta
CurrentModule = SpectraUtils
```

# SpectraUtils

SpectraUtils.jl builds physical spectroscopic models (species states,
transitions, population distributions, and full emission/radiance spectra)
on top of the line-shape/line-synthesis primitives (`LineShape`, `Line`,
`calculate_lines`) provided by the sibling package
[SpectralLines.jl](https://jqfeld.github.io/SpectralLines.jl), which
`SpectraUtils` re-exports.

```@contents
Pages = ["spectra.md", "linewidths.md", "hitran.md"]
Depth = 2
```

## Reference

```@index
```

```@autodocs
Modules = [SpectraUtils]
# Documented on the "Spectra" page; listing them here too duplicates them.
Filter = t -> !(t in (SpectraUtils.SpeciesModel, SpectraUtils.Slab,
                      SpectraUtils.StaticLineshapes, SpectraUtils.EmissionSpectrum,
                      SpectraUtils.RadianceSpectrum, SpectraUtils.AbsorptionCoefficientSpectrum,
                      SpectraUtils.AbsorbanceSpectrum, SpectraUtils.TransmittanceSpectrum))
```
