# Spectra

Five models build a spectrum from a list of quantum
transitions, a population distribution, and a line-shape function:

- [`EmissionSpectrum`](@ref SpectraUtils.EmissionSpectrum): optically-thin
  emission, summing `A_ul · f_u · φ(ν)` per line (photon-rate by default,
  energy-weighted via `amplitude = SpectralRadiance()`).
- [`RadianceSpectrum`](@ref SpectraUtils.RadianceSpectrum): self-absorbing
  radiative transfer through one or more uniform slabs, with an LTE fast path
  (Kirchhoff's law) and a general non-LTE α/j path.
- [`AbsorptionCoefficientSpectrum`](@ref SpectraUtils.AbsorptionCoefficientSpectrum):
  the absorption coefficient `α(ν)` itself.
- [`AbsorbanceSpectrum`](@ref SpectraUtils.AbsorbanceSpectrum): Napierian
  absorbance `α·L = −ln(I/I₀)`.
- [`TransmittanceSpectrum`](@ref SpectraUtils.TransmittanceSpectrum):
  normalized transmission `e^{−αL}`.

Every model takes a `Tuple` of [`SpeciesModel`](@ref SpectraUtils.SpeciesModel)
via its `species` keyword: several species, each with its own line list,
distribution, number density and synthesis method, share one slab, with
coefficients added before anything nonlinear. The flat per-species keywords are
shorthand for a single species. All models but `AbsorptionCoefficientSpectrum`
also take a `slabs` keyword, a `Tuple` of [`Slab`](@ref SpectraUtils.Slab), for a
line of sight through several regions.

## Example

```julia
using SpectraUtils, ComponentArrays

# A toy band: rotational levels J of an excited state, all radiating to one
# ground state. Energies in cm⁻¹.
ground = State(:X, 0, 0, 0.0, 1)
upper(J) = State(:A, 0, J, 20000.0 + 2.0 * J * (J + 1), 2J + 1)
transitions = [ElectricDipoleTransition(ground, upper(J), 1.0e6) for J in 0:20]

# Thermal population of the A state, measured from its lowest level.
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

params = ComponentVector(T = 1500.0, mass = 4.65e-26, gamma = 0.05, sigma_inst = 0.2)
spectrum = sp(wavenumbers, params)
```

Read parameters as `p.T` or `p[Val(:T)]`. Indexing a `ComponentVector` with a
runtime `Symbol` (`p[:T]`) dispatches dynamically and allocates on every line.

## API

```@docs
SpectraUtils.SpeciesModel
SpectraUtils.Slab
SpectraUtils.StaticLineshapes
SpectraUtils.EmissionSpectrum
SpectraUtils.RadianceSpectrum
SpectraUtils.AbsorptionCoefficientSpectrum
SpectraUtils.AbsorbanceSpectrum
SpectraUtils.TransmittanceSpectrum
```

## Synthesis method: binning and wing cutoff

`SpeciesModel`'s `method` field controls how the spectrum is evaluated on the
wavenumber grid: `PointSample`, `BinAverage`, `LDM` or `SparseLDM`, all from
[SpectralLines.jl](https://jqfeld.github.io/SpectralLines.jl).

With `BinAverage(bin_width=...)` each grid point takes the profile's average over a
bin of that width centred on the point, `(1/bin_width) ∫ φ(x) dx`. Set
`bin_width` to the grid spacing. The bins then tile the axis and every line
keeps its area. Point sampling at a spacing comparable to the linewidth misses
most of a narrow line or overweights it, depending on where its centre falls
between grid points. `BinAverage` returns the same quantity in the same units as
`PointSample`, and the two converge as `bin_width` falls below the linewidth. Use
`BinAverage` for `bin_width ≥ 0.3 × FWHM` and `PointSample` (the default) otherwise.

`wing_cut`, a field on both `PointSample` and `BinAverage`, evaluates each line only
within `±wing_cut × hwhm` of its centre, found by binary search on the grid. The
cost per line falls from the whole grid to that window, so the speedup grows
with the grid span over the linewidth. Beyond the cut the line is dropped: a
Lorentzian keeps about `2/(π · wing_cut)` of its area outside `±wing_cut × hwhm`,
0.13% at `wing_cut = 500`.

## Backgrounds

Every model's `background` field takes `nothing`, a constant, or any callable
`(xs, params) -> values`. It receives `params`, so its coefficients are fit
parameters like any other.

[`ChebyshevBackground`](@ref) packages the common case, a polynomial continuum
over a known interval. With [FittingUtils.jl](https://github.com/jqfeld/FittingUtils.jl):

```julia
sp = EmissionSpectrum(; transitions, distribution_function, lineshape_function,
                        background = ChebyshevBackground(306.0, 320.0, 3))

priors = ComponentVector(Trot = Normal(500.0, 300.0),
                         background = [Normal(0.0, 1.0) for _ in 1:4])
res = fit_model(sp, xs, data, priors, LevenbergMarquardtFit())
```

Degree 0 is a constant offset, degree 1 a straight line.

Too low a degree biases the physical parameters, since the fit absorbs the
remaining continuum into the temperature. Take the lowest degree that leaves
no structure in the residuals.
