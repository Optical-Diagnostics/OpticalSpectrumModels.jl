# Linewidth Estimates

Doppler and collisional broadening estimates, returning the width parameters the
line-shape constructors take directly.

## API

```@docs
sigma_doppler
gamma_hard_sphere
gamma_lifetime
gamma_pressure
```

## Example

N₂ at 2330 cm⁻¹ and 300 K, pressure-broadened at 1 atm (all inputs SI except
the wavenumber):

```julia
using SpectraUtils

σ = sigma_doppler(2330.0, 300; m = 28 * 1.66054e-27)          # cm⁻¹, ν̃ in cm⁻¹, m in kg
γ = gamma_hard_sphere(101325, 300;                            # p in Pa
                      μ = 14 * 1.66054e-27,                   # reduced mass, kg
                      cs = 4.5e-19)                           # cross-section, m²
Voigt(σ, γ)
```
