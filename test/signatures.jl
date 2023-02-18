using SparseArrays: SparseVector
using StaticArrays

import GeometricAlgebra:
	dimension,
	basis_vector_norm,
	componentstype

@testset "basis, @basis" begin
	@test length(basis(10)) == 10
	@test length(basis(4, 0:2:4)) == 8
	@test grade(prod(basis(5))) == 5

	@test basis((1,1,1)) == basis("+++")

	let
		@basis "-+++"
		@test v1^2 == v234^2 == -1
	end

	let
		@basis "-+++" allperms=true
		@test v123 == v231 == v312
	end

	let
		@basis (x=1, t=-1) allperms=true
		@test x^2 == xt^2 == tx^2 == -t^2 == 1
	end

	@test all(basis(Multivector{3,1}) .== Multivector.(basis(3)))
	@test all(basis(Multivector{2,0:2}) .== Multivector.(basis(2, 0:2)))
	@test all(eltype.(basis(Multivector{3,0:3,Vector{Bool}})) .=== Bool)

end

@testset "cayleytable" begin
	@test sprint(cayleytable, 1:5, +) == """
	 (↓) + (→) │ 1  2  3  4   5
	───────────┼────────────────
	         1 │ 2  3  4  5   6
	         2 │ 3  4  5  6   7
	         3 │ 4  5  6  7   8
	         4 │ 5  6  7  8   9
	         5 │ 6  7  8  9  10
	"""

	@test sprint(cayleytable, "+++", ∧) == """
	 (↓) ∧ (→) │    1 │   v1     v2    v3 │  v12    v13   v23 │ v123
	───────────┼──────┼───────────────────┼───────────────────┼──────
	         1 │    1 │   v1     v2    v3 │  v12    v13   v23 │ v123
	───────────┼──────┼───────────────────┼───────────────────┼──────
	        v1 │   v1 │    0    v12   v13 │    0      0  v123 │    0
	        v2 │   v2 │ -v12      0   v23 │    0  -v123     0 │    0
	        v3 │   v3 │ -v13   -v23     0 │ v123      0     0 │    0
	───────────┼──────┼───────────────────┼───────────────────┼──────
	       v12 │  v12 │    0      0  v123 │    0      0     0 │    0
	       v13 │  v13 │    0  -v123     0 │    0      0     0 │    0
	       v23 │  v23 │ v123      0     0 │    0      0     0 │    0
	───────────┼──────┼───────────────────┼───────────────────┼──────
	      v123 │ v123 │    0      0     0 │    0      0     0 │    0
	"""
end

# Test that things generally work with Multivectors with various other component array types
# Fixed size array types such as MVector or SVector may be most performant in few dimensions,
# while something like SparseVector is suited to very large algebras.
struct SigWithCompsType{Sig,S} end
dimension(::SigWithCompsType{Sig}) where {Sig} = dimension(Sig)
basis_vector_norm(::SigWithCompsType{Sig}, i) where {Sig} = basis_vector_norm(Sig, i)

componentstype(::SigWithCompsType{Sig,<:Vector}, N, T) where {Sig} = Vector{T}
componentstype(::SigWithCompsType{Sig,<:MVector}, N, T) where {Sig} = MVector{N,T}
componentstype(::SigWithCompsType{Sig,<:SVector}, N, T) where {Sig} = SVector{N,T}
componentstype(::SigWithCompsType{Sig,<:SparseVector}, N, T) where {Sig} = SparseVector{T}

# SparseVectors are appropriate for large dims, but symbolic optimization only for small dim — you wouldn’t have both
GeometricAlgebra.use_symbolic_optim(::SigWithCompsType{Sig,<:SparseVector}) where {Sig} = false

@testset "other component array types" begin
	for C in [Vector, MVector, SVector, SparseVector]

		sig = SigWithCompsType{3,C}()
		b = basis(sig)

		@test b[1] + b[2] isa Multivector{sig,1,<:C}
		@test exp(b[1]b[2]) isa Multivector{sig,0:2:2,<:C}

		m = 1 + b[1]
		@test m isa Multivector{sig,K,<:C} where K
		@test grade(m, 1) + b[2] isa Multivector{sig,1,<:C}

		u = sum(b)
		@test u/u isa Multivector{sig,K,<:C} where K
		@test sin(~u) isa Multivector{sig,K,<:C} where K
	end
end