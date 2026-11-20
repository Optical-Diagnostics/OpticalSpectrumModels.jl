"""
    make_transition_list(source) -> Vector{ElectricDipoleTransition}

The transitions of a line-list `source`, the input of every spectrum model.
Methods exist for a `Vector{HITRAN.Line}` and, with `using ExoMol`, for an
`ExoMol.Isotopologue`. Add a method to support another source.
"""
make_transition_list(x) = throw(ArgumentError(
    "make_transition_list has no method for $(typeof(x)). It accepts a " *
    "Vector{HITRAN.Line}, or an ExoMol.Isotopologue after `using ExoMol`."))
