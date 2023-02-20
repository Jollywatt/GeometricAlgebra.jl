using GeometricAlgebra:
	isscalar,
	componentbits,
	componentindex,
	componentindices
using GeometricAlgebra.StaticArrays
using GeometricAlgebra.SparseArrays

@testset "constructors" begin
	@test zero(BasisBlade{2}(42, 0b11)) === BasisBlade{2,0,Int}(0, 0b00)
	@test one(BasisBlade{(1,1)}(42.0, 0b11)) === BasisBlade{(1,1),0,Float64}(1, 0b00)

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

	@test zero(Multivector{3,1}) isa Multivector
	@test zero(Multivector{3,1}, Float32) |> eltype === Float32
	@test zero(Multivector{3,2}, Float32) |> grade === 2
end

@testset "components" begin
	for n in 0:6, k in [0:n..., 0:n, 0:2:n, n > 0 ? (0, n) : 0]
		a = Multivector{n,k}
		@test componentindex.(a, componentbits(a)) == 1:ncomponents(a)

		if k isa Integer
			bits = componentbits(a)[componentindices(a, k)]
			@test length(bits) == binomial(n, k)
			@test all(count_ones.(bits) .== k)
		end
	end
end

@testset "iseven/isodd" begin
	@test iseven(BasisBlade{3}(42, 0b11))
	@test isodd(BasisBlade{3}(42, 0b111))
	@test iseven(Multivector{4,0:2:4}(ones(8)))
	@test isodd(Multivector{4,1:2:4}(ones(8)))

	m = Multivector{3,0:3}(ones(8))
	@test !iseven(m) && !isodd(m)
end

@testset "scalar" begin
	@test scalar(BasisBlade{3}(42)) == 42
	@test scalar(BasisBlade{3}(42, 0b111)) == 0
	@test scalar(Multivector{3,0:3}(1:8)) == 1
	@test scalar(Multivector{3,3}([42])) == 0
end
