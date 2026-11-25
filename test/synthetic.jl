# Two synthetic transitions on distinct vibrational levels, shared by the
# model tests. `bd` is a Boltzmann distribution over exactly these states.
using OpticalSpectrumModels: h, c0, c0_cm, kB

const c2 = 100 * h * c0 / kB   # cm·K, second radiation constant

const lower0 = State(:X, 0, nothing, 0.0, 1.0)
const upper0 = State(:X, 0, nothing, 1000.0, 3.0)
const t1     = ElectricDipoleTransition(lower0, upper0, 1000.0, 10.0)

const lower1 = State(:X, 1, nothing, 0.0, 1.0)
const upper1 = State(:X, 1, nothing, 1050.0, 5.0)
const t2     = ElectricDipoleTransition(lower1, upper1, 1050.0, 8.0)

const transitions = [t1, t2]

Qfun(T) = 1.0 + 3.0 * exp(-c2 * 1000.0 / T) + 5.0 * exp(-c2 * 1050.0 / T)
const bd = Boltzmann(Qfun, 0.0, State(:X))
