```@meta
CurrentModule = OpticalSpectrumModels
```

# OpticalSpectrumModels

OpticalSpectrumModels.jl builds physical spectroscopic models (species states,
transitions, population distributions, and full emission/radiance spectra)
on top of the line-shape/line-synthesis primitives (`LineShape`, `Line`,
`calculate_lines`) provided by the sibling package
[SpectralLines.jl](https://jqfeld.github.io/SpectralLines.jl), which
`OpticalSpectrumModels` re-exports.

```@contents
Pages = ["spectra.md", "linewidths.md", "hitran.md"]
Depth = 2
```

## Reference

```@index
```

```@autodocs
Modules = [OpticalSpectrumModels]
# Documented on the "Spectra" page; listing them here too duplicates them.
Filter = t -> !(t in (OpticalSpectrumModels.SpeciesModel, OpticalSpectrumModels.Slab,
                      OpticalSpectrumModels.StaticLineshapes, OpticalSpectrumModels.EmissionSpectrum,
                      OpticalSpectrumModels.RadianceSpectrum, OpticalSpectrumModels.AbsorptionCoefficientSpectrum,
                      OpticalSpectrumModels.AbsorbanceSpectrum, OpticalSpectrumModels.TransmittanceSpectrum))
```
