#= GeometricAlgebra.jl
An implementation of geometric algebras (a.k.a. Clifford algebras) in Julia.
=#

module GeometricAlgebra

import Base: ==, *, /, \, +, -, ^, ~
import SparseArrays: SparseVector, spzeros

export AbstractMultivector
export Blade, Multivector, MixedMultivector

export signature, dim, basis, vol
export grade, scalar, isscalar
export blades

export reversion, ~
export hodgedual, ★

export dot, ⋅
export wedge, ∧
export scalar_prod, ∗
export contractr, ⨽
export contractl, ⨼

export EuclideanSignature, OffsetSignature

include("metric-signature.jl")

include("ublade.jl")

include("multivector.jl")

include("algebra.jl")

include("show.jl")

end # module
