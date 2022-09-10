"""
Multivectors
---

Implementation of multivectors, or ``k``-vectors, as in geometric (Clifford) algebra.
"""
module Multivectors

using StaticArrays

export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector
export dimension, grade
export geometric_prod, scalar_prod, wedge, ∧, reversion

include("bits.jl")
include("types.jl")
include("utils.jl")
include("signatures.jl")
include("algebra.jl")
include("special.jl")
include("show.jl")


end # module Multivectors
