"""
    OpticalSpectrumModels

Utilities for building physical spectroscopic models (species states,
transitions, population distributions, and full emission/radiance spectra)
on top of the line-shape/line-synthesis primitives provided by the sibling
package `SpectralLines.jl`.
"""
module OpticalSpectrumModels

using SpectralLines
import SpectralLines: convolve, line_index_range, spectrum_eltype, LineShape, numtype_convert, check_axes
export fwhm, hwhm, integrate, convolution_support
export Lorentzian, Gaussian, Voigt, PseudoVoigt, Triangular
export PseudoVoigtModel, Thompson, Kielkopf, Olivero, LiuLin, Rocco
export Line, update_lines!, calculate_lines, calculate_lines!, LineSynthesisMethod
export PointSample, BinAverage, LDM, SparseLDM, LDMWorkspace, SparseLDMWorkspace
export convolve


include("./constants.jl")
include("./import.jl")

export HITRAN
include("./HITRAN/HITRAN.jl")

export sigma_doppler, gamma_hard_sphere, gamma_lifetime, gamma_pressure
export molecular_mass_kg
include("./linewidths.jl")

export convolve_spectrum
include("./convolution.jl")

export LineCache, fill_lines!
include("./line_cache.jl")

export electronic_state, vibrational_state, vibronic_state, rotational_state, degeneracy, energy, State
export ElectricDipoleTransition, make_transition_list
export upper_manifold
include("./transitions.jl")

include("./hitran_interop.jl")

export DistributionFunction, StateResolved, Boltzmann, ProductDistribution, ConstantDistribution, VariableDistribution, is_LTE
include("./distribution_functions.jl")

export ChebyshevBackground, n_coefficients, coefficient_key
include("./backgrounds.jl")

export EmissionSpectrum, RadianceSpectrum, planck_radiance
export AbsorptionCoefficientSpectrum, AbsorbanceSpectrum, TransmittanceSpectrum
export SpeciesModel, Slab, StaticLineshapes
export PhotonRate, SpectralRadiance
export calculate_radiance_lines, calculate_radiance_lines!
export emission_line, absorption_line, radiance_line
include("spectrum_models.jl")


end
