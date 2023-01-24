using SparseArrays: SparseVector
using StaticArrays: StaticVector

@testset "basis, @basis" begin
	@test length(basis(10)) == 10
	@test length(basis(4, grade=0:2:4)) == 8
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

struct SigWithStorageType{Sig,S} end
GeometricAlgebra.dimension(::SigWithStorageType{Sig}) where {Sig} = dimension(Sig)
GeometricAlgebra.basis_vector_norm(::SigWithStorageType{Sig}, i) where {Sig} = GeometricAlgebra.basis_vector_norm(Sig, i)
GeometricAlgebra.componentstype(::SigWithStorageType{Sig,<:Vector}, N, T) where {Sig} = Vector{T}
GeometricAlgebra.componentstype(::SigWithStorageType{Sig,<:MVector}, N, T) where {Sig} = isbitstype(T) ? MVector{N,T} : Vector{T} # MVectors only support setindex! for isbits types
GeometricAlgebra.componentstype(::SigWithStorageType{Sig,<:SparseVector}, N, T) where {Sig} = SparseVector{T}


@testset "other storage types" begin
	for C in [Vector, SparseVector, MVector]

		sig = SigWithStorageType{3,C}()
		b = basis(sig)

		@test b[1] + b[2] isa Multivector{sig,1,<:C}
		@test exp(b[1]b[2]) isa Multivector{sig,0:2:2,<:C}

		m = 1 + b[1]
		@test m isa Multivector{sig,K,<:C} where K
		@test grade(m, 1) + b[2] isa Multivector{sig,1,<:C}
	end
end