"""
    sigma_doppler(f0, T; m) -> σ [cm⁻¹]

Gaussian standard deviation of the Doppler-broadened line at wavenumber `f0`
[cm⁻¹], temperature `T` [K], and molecular mass `m` [kg]: the parameter
`Gaussian`/`Voigt` expect (FWHM is `√(8 ln 2) · σ`, the 1/e half-width `√2 · σ`).

```math
\\sigma = \\frac{f_0}{c}\\sqrt{\\frac{k_B T}{m}}
```
"""
@inline sigma_doppler(f0, T; m) = sqrt(kB) / c0 * sqrt(T / m) * f0

"""
    sigma_doppler(f0, T, species) -> σ [cm⁻¹]

As above, with the mass from `species` instead of a supplied `m`. Uses the true
isotopologue mass, not the sum of mass numbers.

Works for any `species` with a `molecular_mass_kg` method:
`ExoMol.Isotopologue` via `SpectraUtilsExoMolExt`, or for HITRAN
`molecular_mass_kg(HITRAN.molparam(mol_id, iso))`.
"""
@inline sigma_doppler(f0, T, species) = sigma_doppler(f0, T; m = molecular_mass_kg(species))

"""
    molecular_mass_kg(species) -> m [kg]

Molecular mass of `species` in kilograms. Implemented for HITRAN `molparam`
entries here and for `ExoMol.Isotopologue` in `SpectraUtilsExoMolExt`. Extend it
for other line-list sources to make [`sigma_doppler`](@ref)'s three-argument
form work with them.
"""
function molecular_mass_kg end

# HITRAN's molparam table gives a molar mass in g/mol.
@inline molecular_mass_kg(e::@NamedTuple{isoname::String, abundance::Float64,
                                         Q296K::Float64, gj::Int, molar_mass::Float64}) =
    e.molar_mass * 1.0e-3 / N_A

"""
    gamma_hard_sphere(p, T; μ, cs) -> γ [cm⁻¹]

Lorentzian HWHM of pressure broadening at pressure `p` [Pa], temperature `T` [K],
reduced mass `μ` [kg], and collision cross-section `cs` [m²], using the hard-sphere
model: the collision rate `n·cs·v̄` over `2πc`,

```math
\\gamma = \\frac{n\\,\\sigma\\,\\bar{v}}{2\\pi c},\\qquad
n = \\frac{p}{k_B T},\\quad \\bar{v} = \\sqrt{\\frac{8 k_B T}{\\pi \\mu}}
```

N₂ self-broadening at 1 atm and 300 K (`cs ≈ 4.5e-19` m²) gives γ ≈ 0.04 cm⁻¹.
"""
@inline gamma_hard_sphere(p, T; μ, cs) = p / kB / T * cs * sqrt(8 / pi * kB * T / μ) / (2π * c0 * 100)

"""
    gamma_lifetime(tau) -> γ [cm⁻¹]

Lorentzian HWHM due to a finite state lifetime `tau` [s]:

```math
\\gamma = \\frac{1}{4\\pi c\\,\\tau}
```

where ``c`` is the speed of light in cm s⁻¹. For a purely radiative state,
``\\tau = 1 / \\sum_l A_{ul}``; for predissociating states it should include
the non-radiative decay rate as well.

Add to the pressure-broadening ``\\gamma`` before constructing a `Voigt` profile:
```julia
γ_total = gamma_hard_sphere(p, T; μ, cs) + gamma_lifetime(tau)
Voigt(sigma_doppler(f0, T; m), γ_total)
```
"""
@inline gamma_lifetime(tau) = inv(4 * pi * c0_cm * tau)

"""
    gamma_pressure(p_atm, T; gamma_ref, n, Tref=296.0) -> γ [cm⁻¹]

Lorentzian HWHM of pressure broadening from a tabulated HITRAN/ExoMol
coefficient `gamma_ref` [cm⁻¹/atm] at reference temperature `Tref` [K] (296 K by
convention), scaled to pressure `p_atm` [atm] and temperature `T` [K]:

```math
\\gamma = \\gamma_{ref}\\,p_{atm}\\,(T_{ref}/T)^n
```

The tabulated-coefficient approach of HITRAN/ExoMol line lists, as opposed to
[`gamma_hard_sphere`](@ref)'s kinetic-theory model. `gamma_ref`/`n` come from a
per-line `.broad` file, or from a dataset-wide default in the `.def.json` for
datasets without one (SiO EBJT: `gamma_ref=0.07`, `n=0.5`).

For a mixture, sum one term per diluent weighted by partial pressure:

```julia
γ = gamma_pressure(p_atm * (1 - x_self), T; gamma_ref=γ_air, n=n_air) +
    gamma_pressure(p_atm * x_self, T; gamma_ref=γ_self, n=n_self)
```
"""
@inline gamma_pressure(p_atm, T; gamma_ref, n, Tref=296.0) = gamma_ref * p_atm * (Tref / T)^n

"""
    gamma_pressure(p_atm, T, x_self; gamma_air, n_air,
                            gamma_self=gamma_air, n_self=n_air, Tref=296.0) -> γ [cm⁻¹]

Air plus self broadening in one call: the two-term sum above, with partial
pressures derived from the mole fraction `x_self`.

`gamma_self`/`n_self` default to the air values, matching ExoMol datasets that
declare one dataset-wide coefficient rather than separate air and self ones
(SiO EBJT).

```julia
γ = gamma_pressure(P / 101325, T, x; gamma_air = 0.07, n_air = 0.5)
```
"""
@inline function gamma_pressure(p_atm, T, x_self;
                                         gamma_air, n_air,
                                         gamma_self = gamma_air, n_self = n_air,
                                         Tref = 296.0)
    return gamma_pressure(p_atm * (1 - x_self), T; gamma_ref = gamma_air, n = n_air, Tref) +
           gamma_pressure(p_atm * x_self, T; gamma_ref = gamma_self, n = n_self, Tref)
end

