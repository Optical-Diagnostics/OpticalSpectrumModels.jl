# README

Data used for comparisons with calculations from `OpticalSpectrumModels.jl`

## Radis
Radis version: 0.17
SiO example from the Radis documentation:

```python 
from radis import calc_spectrum

s = calc_spectrum(
    1080,
    1320,  # cm-1
    species="SiO",
    isotope="1",
    pressure=1.01325,  # bar
    Tgas=1000,  # K
    mole_fraction=0.1,
    path_length=1,  # cm
    databank=(
        "exomol",
        "EBJT",
    ),  # Or simply use 'exomol' to auto-select the recommended database (if one exists)
)
s.apply_slit(1, "cm-1")  # simulate an experimental slit

s.savetxt("radis_SiO_example_radiance.csv", "radiance")
s.savetxt("radis_SiO_example_absorbance.csv", "absorbance")

```

## OES toolbox
Instrument profile (pseudo-Voigt):
- Shape (0-1): 0.5
- Width (nm): 0.01

Transition:
- N2(C-B)
- Trot: 300 K
- Tvib: 1200 K

