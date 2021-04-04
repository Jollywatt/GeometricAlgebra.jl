module GeometricAlgebra

import Base: ==, *, /, \, +, -, ^
import InteractiveUtils: subtypes
import Combinatorics: parity

# for use in documentation to give module-independent cross references
fullname(method) = "$(parentmodule(method)).$(nameof(method))"


include("unit_blades.jl")


export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector
include("types.jl")

include("signatures.jl")

export ∧, ⋅
export ~, reversion, involute
export grade, grades
include("algebra.jl")

export basis
include("convenience.jl")

include("show.jl")


export best_type

end # module
