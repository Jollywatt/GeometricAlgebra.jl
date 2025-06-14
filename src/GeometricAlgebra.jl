module GeometricAlgebra

using StaticArrays, SparseArrays
using Combinatorics: permutations, parity
using PrettyTables: pretty_table
using Random: AbstractRNG

export SVector, MVector

export AbstractMultivector, BasisBlade, Multivector, Grade
export basis, @basis, cayleytable
export BasisDisplayStyle
export signature, dimension, ncomponents, grade, eachgrade, blades
export Cl
export scalar, isscalar, ishomogeneous, isblade
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
export matrix_repr, outermorphism, embed
export factorblade
export @symbolicga
export MiniCAS

include("MiniCAS/MiniCAS.jl")

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
include("promotion.jl")

include("Conformal/cga.jl")

end # module GeometricAlgebra
