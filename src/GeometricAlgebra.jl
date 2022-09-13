"""
GeometricAlgebra
---

Implementation of multivectors, or ``k``-vectors, as in geometric (Clifford) algebra.
"""
module GeometricAlgebra

using StaticArrays

export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector
export basis
export signature, dimension, grade
export geometric_prod, scalar_prod, wedge, ∧, reversion, involution

include("bits.jl")
include("types.jl")
include("utils.jl")
include("signatures.jl")
include("algebra.jl")
include("special.jl")
include("generated.jl")
include("show.jl")


end # module GeometricAlgebra
