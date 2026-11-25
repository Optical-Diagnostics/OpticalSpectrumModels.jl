using OpticalSpectrumModels
using Test

# CODATA values, written out so a wrong package constant cannot cancel.
const kB = 1.380649e-23    # J/K
const c  = 299792458.0     # m/s

@testset "gamma_hard_sphere" begin
    # N₂ self-broadening at 1 atm, 300 K: γ = n·σ·v̄ / (2πc) ≈ 0.04 cm⁻¹.
    p, T = 101325.0, 300.0
    μ  = 2.33e-26   # kg, reduced mass of N₂–N₂
    cs = 4.5e-19    # m², hard-sphere cross-section
    n = p / (kB * T)
    v̄ = sqrt(8 * kB * T / (pi * μ))
    @test gamma_hard_sphere(p, T; μ, cs) ≈ n * cs * v̄ / (2π * c * 100) rtol = 1e-12
end

@testset "gamma_lifetime" begin
    # γ = 1/(4πcτ) in cm⁻¹.
    @test gamma_lifetime(1e-9) ≈ inv(4π * c * 100 * 1e-9) rtol = 1e-12
end

@testset "sigma_doppler" begin
    # σ = ν̃₀·√(kBT/m)/c.
    m = 28.0 * 1.66054e-27
    @test sigma_doppler(2330.0, 296.0; m) ≈ 2330.0 / c * sqrt(kB * 296.0 / m) rtol = 1e-12

    # The species forms read the mass from a HITRAN molparam entry.
    e = OpticalSpectrumModels.HITRAN.molparam(2, 1)          # CO2, main isotopologue
    @test molecular_mass_kg(e) ≈ e.molar_mass * 1e-3 / 6.02214076e23 rtol = 1e-12
    @test sigma_doppler(1000.0, 300.0, e) == sigma_doppler(1000.0, 300.0; m = molecular_mass_kg(e))
end

@testset "gamma_pressure mixture overload" begin
    p_atm, T, x = 1.0, 1000.0, 0.1
    γa, na, γs, ns = 0.07, 0.5, 0.09, 0.6

    manual = OpticalSpectrumModels.gamma_pressure(p_atm * (1 - x), T; gamma_ref = γa, n = na) +
             OpticalSpectrumModels.gamma_pressure(p_atm * x,       T; gamma_ref = γs, n = ns)
    @test OpticalSpectrumModels.gamma_pressure(p_atm, T, x;
              gamma_air = γa, n_air = na, gamma_self = γs, n_self = ns) ≈ manual rtol = 1e-14

    # Self defaults to air (the single-coefficient ExoMol case), and then the
    # mixture reduces to the pure-pressure result regardless of x.
    both = OpticalSpectrumModels.gamma_pressure(p_atm, T, x; gamma_air = γa, n_air = na)
    @test both ≈ OpticalSpectrumModels.gamma_pressure(p_atm, T; gamma_ref = γa, n = na) rtol = 1e-14
end
