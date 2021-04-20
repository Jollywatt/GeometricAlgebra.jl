module GeometricAlgebra

import Base: ==, *, /, \, +, -, ^
import InteractiveUtils: subtypes
import Combinatorics: parity, powerset, permutations
import OffsetArrays: OffsetArray # only ever used by basis(::OffsetSignature)
using SparseArrays
using StaticArrays
using TaylorSeries: Taylor1


include("unit_blades.jl")

export AbstractMultivector, HomogeneousMultivector
export Blade, Multivector, MixedMultivector
include("types.jl")

export signature, dimension
export OffsetSignature
include("signatures.jl")

export ∧, ⋅, ⨼, ⨽, ∗
export ~, reversion, involute
export grade, grades, isscalar, scalar
export unit_pseudoscalar, dual
include("algebra.jl")

include("special.jl")

export basis, @basis, @basisall
include("convenience.jl")

include("show.jl")


export best_type

end # module
