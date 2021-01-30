#= GeometricAlgebra.jl
An implementation of geometric algebras (a.k.a. Clifford algebras) in Julia.
=#

module GeometricAlgebra

import Base: ==, *, /, \, +, -, ^, ~
import SparseArrays: SparseVector, spzeros
import Combinatorics: powerset, permutations

export AbstractMultivector
export Blade, Multivector, MixedMultivector

export signature, dim, basis, vol
export grade, grades
export scalar, isscalar
export blades

export @basis

export reversion, ~

export dot, ⋅
export wedge, ∧
export contractr, ⨽
export contractl, ⨼

export EuclideanSignature, OffsetSignature

include("metric-signature.jl")

include("ublade.jl")

include("multivector.jl")

include("algebra.jl")

include("show.jl")

include("convenience.jl")

end # module
