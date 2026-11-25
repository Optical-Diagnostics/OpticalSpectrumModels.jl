using OpticalSpectrumModels
using ExoMol
using ComponentArrays

# N₂ first positive system: B³Πg → A³Σu⁺
# Spectral range: ~10000–19500 cm⁻¹ (~510–1000 nm)
const _N2 = load_isotopologue("N2", "14N2", "WCCRMT"; wn_range=(10000.0, 19500.0))

tl = make_transition_list(_N2)
@debug "Total transitions in range" length(tl)

# The electronic state labels of this dataset
@info "Electronic states (upper)" collect(unique(electronic_state(t.upper) for t in tl))
@info "Electronic states (lower)" collect(unique(electronic_state(t.lower) for t in tl))

# First positive system: B³Πg (upper) → A³Σu⁺ (lower)
# State labels from WCCRMT: :B3Pig for B³Πg, Symbol("A3Sigmau+") for A³Σu⁺.
filter!(t -> electronic_state(t.upper) == Symbol("B3Pig"), tl)
filter!(t -> electronic_state(t.lower) == Symbol("A3Sigmau+"), tl)
@info "First positive system transitions" length(tl)

# Boltzmann population over all B³Πg vibrational levels
const bd = Boltzmann(_N2, State(Symbol("B3Pig")); Ts=200:100:6000)


sp = EmissionSpectrum(;
    species = SpeciesModel(;
        transitions           = tl,
        distribution_function = bd,
        lineshape_function    = (t, p) -> Voigt(sigma_doppler(t.wavenumber, p.T, _N2), 0.05),
        method                = PointSample(wing_cut=100.0)),
)

xs     = collect(range(10000.0, 19500.0; length=2000))
params = ComponentVector(T=5000.0)

spectrum = sp(xs, params)
@info "Spectrum computed" length=length(xs) peak=maximum(spectrum)
