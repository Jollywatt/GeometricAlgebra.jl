"""
Multivectors
---

Implementation of multivectors, or ``k``-vectors, as in geometric (Clifford) algebra.
"""
module Multivectors

import Base: ==

export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector

include("bits.jl")
include("types.jl")
include("signatures.jl")


end # module Multivectors
