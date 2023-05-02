module GeometricAlgebra

using StaticArrays, SparseArrays
using Combinatorics: permutations, parity
using SymbolicUtils
using PrettyTables: pretty_table

export MVector

export AbstractMultivector, BasisBlade, Multivector
export basis, @basis, cayleytable
export signature, dimension, ncomponents, grade, eachgrade, blades
export Cl, @Cl_str
export scalar, isscalar, ishomogeneous
export geometric_prod,
	scalar_prod, ⊙,
	wedge, ∧,
	antiwedge, ∨,
	inner, ⋅,
	lcontract, ⨼,
	rcontract, ⨽,
	clifford_conj, var"'ᶜ",
	sandwich_prod
export reversion, involution
export flipdual, hodgedual, invhodgedual, ldual, rdual
export matrix_repr, outermorphism

include("bits.jl")
include("types.jl")
include("grades.jl")
include("signatures.jl")
include("basis.jl")
include("generated.jl")
include("algebra.jl")
include("special.jl")
include("show.jl")
include("utils.jl")

end # module GeometricAlgebra
