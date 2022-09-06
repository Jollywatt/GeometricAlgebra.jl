"""
Multivectors
---

Implementation of multivectors, or ``k``-vectors, as in geometric (Clifford) algebra.
"""
module Multivectors

export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector
export dimension, grade


include("bits.jl")
include("types.jl")
include("signatures.jl")
include("algebra.jl")
include("utils.jl")


end # module Multivectors
