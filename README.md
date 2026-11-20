# SpectraUtils

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jqfeld.github.io/SpectraUtils.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jqfeld.github.io/SpectraUtils.jl/dev/)
[![Build Status](https://github.com/jqfeld/SpectraUtils.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jqfeld/SpectraUtils.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jqfeld/SpectraUtils.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jqfeld/SpectraUtils.jl)

SpectraUtils.jl builds physical spectroscopic models (species states,
transitions, population distributions, and full emission/radiance spectra) on
top of the line-shape and line-synthesis primitives of
[SpectralLines.jl](https://github.com/jqfeld/SpectralLines.jl), which it
re-exports. It also reads HITRAN/HITEMP line lists and estimates Doppler and
collisional linewidths from physical parameters.

## Example

A toy band of rotational lines at 1500 K, Doppler- and pressure-broadened, seen
through a Gaussian instrument:

```julia
using SpectraUtils, ComponentArrays

ground = State(:X, 0, 0, 0.0, 1)
upper(J) = State(:A, 0, J, 20000.0 + 2.0 * J * (J + 1), 2J + 1)
transitions = [ElectricDipoleTransition(ground, upper(J), 1.0e6) for J in 0:20]

c2 = 1.4387769                                   # hc/k, cm·K
Q(T) = sum((2J + 1) * exp(-c2 * 2.0 * J * (J + 1) / T) for J in 0:20)
bd = Boltzmann(Q, 20000.0, State(:A))

wavenumbers = range(19990.0, 20900.0, length = 9101)
sp = EmissionSpectrum(;
    species = SpeciesModel(;
        transitions,
        distribution_function = bd,
        lineshape_function = (t, p) -> Voigt(sigma_doppler(t.wavenumber, p.T; m = p.mass), p.gamma),
        method = BinAverage(bin_width = step(wavenumbers), wing_cut = 50.0)),
    instrument_function = (wn, p) -> Gaussian(p.sigma_inst),
)
spectrum = sp(wavenumbers, ComponentVector(T = 1500.0, mass = 4.65e-26, gamma = 0.05, sigma_inst = 0.2))
```

See the [documentation](https://jqfeld.github.io/SpectraUtils.jl/dev/) for
self-absorbing, absorption and multi-slab models, HITRAN and ExoMol line lists,
and backgrounds.
