module GeometricAlgebra

using StaticArrays, SparseArrays
using Combinatorics: permutations, powerset
using SymbolicUtils
using PrettyTables: pretty_table

export AbstractMultivector, BasisBlade, Multivector
export basis, @basis, @basisall, cayleytable
export signature, dimension, ncomponents, grade, eachgrade, blades
export scalar, isscalar, ishomogeneous
export geometric_prod, scalar_prod, ⊙, wedge, ∧, antiwedge, ∨, inner, ⋅, lcontract, ⨼, rcontract, ⨽, clifford_conj, var"'ᶜ"
export reversion, involution
export flipdual, hodgedual, poincaredual
export matrix_repr
export Cl

include("bits.jl")
include("types.jl")
include("grades.jl")
include("signatures.jl")
include("generated.jl")
include("algebra.jl")
include("special.jl")
include("show.jl")
include("utils.jl")

end # module GeometricAlgebra
