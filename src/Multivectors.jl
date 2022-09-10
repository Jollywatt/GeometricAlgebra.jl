"""
Multivectors
---

Implementation of multivectors, or ``k``-vectors, as in geometric (Clifford) algebra.
"""
module Multivectors

export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector
export dimension, grade
export geometric_prod, scalar_prod, wedge, ∧

include("bits.jl")
include("types.jl")
include("utils.jl")
include("signatures.jl")
include("algebra.jl")
include("show.jl")


end # module Multivectors
