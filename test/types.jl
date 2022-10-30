using GeometricAlgebra:
	isscalar,
	componentbits,
	componentindex,
	componentslice,
	unify_grades,
	resulting_multivector_type
using GeometricAlgebra.StaticArrays
using GeometricAlgebra.SparseArrays

@testset "constructors" begin
	@test zero(BasisBlade{2}(0b11 => 42)) === BasisBlade{2,0,Int}(0b00, 0)
	@test one(BasisBlade{(1,1)}(0b11 => 42.0)) === BasisBlade{(1,1),0,Float64}(0b00, 1)

	for T in [
		Int,
		BasisBlade{2,0,Bool},
		Multivector{3,1,Vector{Int}},
		Multivector{3,0:3,Vector{Float64}},
		Multivector{(-1,+1,+1,+1),0:2:4,Vector{Float64}},
	]
		@test zero(T) isa T
		@test zero(zero(T)) isa T
		@test iszero(zero(T))
		@test isone(one(T))
		@test isscalar(zero(T))
	end
end

@testset "components" begin
	for n in 0:6, k in [0:n..., 0:n, 0:2:n, n > 0 ? (0, n) : 0]
		a = Multivector{n,k}
		@test componentindex.(a, componentbits(a)) == 1:ncomponents(a)

		if k isa Integer
			bits = componentbits(a)[componentslice(a, k)]
			@test length(bits) == binomial(n, k)
			@test all(count_ones.(bits) .== k)
		end
	end
end

@testset "grade inference" begin
	for n in 0:4,
		p in [0:n..., 0:n, 0:2:n, n > 0 ? (0, n) : 0],
		q in [0:n..., 0:n, 0:2:n, n > 0 ? (0, n) : 0]

		pq = unify_grades(n, p, q)
		@test p ⊆ pq ⊇ q
		@test length(pq) != 1 || pq isa Integer

		a = resulting_multivector_type(+, Multivector{n,p}, Multivector{n,q})
		@test p ⊆ grade(a) ⊇ q
	end
end