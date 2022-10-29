module GeometricAlgebra

using StaticArrays, SparseArrays
using Combinatorics: permutations, powerset
using SymbolicUtils
using PrettyTables: pretty_table

export AbstractMultivector, BasisBlade, Multivector
export basis, @basis, @basisall, cayleytable
export signature, dimension, grade, scalar, ncomponents, blades
export geometric_prod, scalar_prod, ⊙, wedge, ∧, inner, ⋅, lcontract, ⨼, rcontract, ⨽, clifford_conj, var"'ᶜ"
export reversion, involution
export flipdual, hodgedual, poincaredual
export vector_repr, matrix_repr
export Cl

include("bits.jl")
include("types.jl")
include("utils.jl")
include("signatures.jl")
include("algebra.jl")
include("special.jl")
include("generated.jl")
include("show.jl")

end # module GeometricAlgebra
