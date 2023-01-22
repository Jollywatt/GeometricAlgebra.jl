using GeometricAlgebra:
	isscalar,
	componentbits,
	componentindex,
	componentslice
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
		@test isscalar(one(T))
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

@testset "iseven/isodd" begin
	@test iseven(BasisBlade{3}(0b11 => 42))
	@test isodd(BasisBlade{3}(0b111 => 42))
	@test iseven(Multivector{4,0:2:4}(ones(8)))
	@test isodd(Multivector{4,1:2:4}(ones(8)))

	m = Multivector{3,0:3}(ones(8))
	@test !iseven(m) && !isodd(m)
end

@testset "scalar" begin
	@test scalar(BasisBlade{3}(0 => 42)) == 42
	@test scalar(BasisBlade{3}(0b111 => 42)) == 0
	@test scalar(Multivector{3,0:3}(1:8)) == 1
	@test scalar(Multivector{3,3}([42])) == 0
end

@testset "similar" begin
	M = similar(Multivector{4,1}, BasisBlade{1,0}(0b0 => one(Float32)))
	@test signature(M) == 4 && grade(M) == 1 && eltype(M) == Float32
end