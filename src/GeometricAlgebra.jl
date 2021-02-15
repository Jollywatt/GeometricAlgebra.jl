#= GeometricAlgebra.jl
An implementation of geometric algebras (a.k.a. Clifford algebras) in Julia.
=#

module GeometricAlgebra

import Base: ==, *, /, \, +, -, ^, ~

using InteractiveUtils: subtypes
using OffsetArrays
using SparseArrays: SparseVector, spzeros
using TaylorSeries: Taylor1
using Combinatorics: powerset, permutations


export AbstractMultivector, Blade, Multivector, MixedMultivector

export basis, @basis, @basisperm
export signature, dimension
export grade, grades
export scalar, isscalar
export vol
export blades

export reversion, ~
export wedge, ∧
export dot, ⋅
export contractr, ⨽
export contractl, ⨼

export EuclideanSignature, OffsetSignature

include("metric-signature.jl")
include("ublade.jl")
include("multivector.jl")
include("show.jl")
include("algebra.jl")
include("convenience.jl")

end # module
