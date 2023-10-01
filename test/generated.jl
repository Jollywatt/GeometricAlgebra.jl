using GeometricAlgebra:
	zeroslike,
	make_symbolic,
	SingletonVector
using GeometricAlgebra.SymbolicUtils


@testset "symbolic components" begin
	@test zeroslike(Vector{Int}, 1) == [0]
	@test zeroslike(Vector{Any}, 2) == [0, 0]
	@test zeroslike(Vector{SymbolicUtils.Symbolic}, 3) == [0, 0, 0]
	@test zeroslike(Vector{Union{Int,SymbolicUtils.Symbolic}}, 4) == [0, 0, 0, 0]


	b = BasisBlade{3}(SymbolicUtils.Sym{Real}(:b), 0b101)
	kv = make_symbolic(Multivector{3,1,Vector{Int}}, :kv)
	mv = make_symbolic(Multivector{3,0:3,Vector{Int}}, :mv)

	for a in [b, kv, mv]
		@test !iszero(a)
		@test !isone(a)

		@test iszero(zero(a))
		@test isone(one(a))

		@test Multivector(a) isa Multivector
	end

	@testset "symbolic algebra" begin
		@test b + 2kv + π*mv isa Multivector
		@test grade(b^2) == 0
		@test kv*5mv isa Multivector
		@test 1 + kv isa Multivector

		@test kv.^(0:3) isa Vector{<:AbstractMultivector}

		@test iszero(kv∧kv)
		@test b∧kv∧mv isa Multivector
	end

end

@testset "SingletonVector" begin
	@test collect(SingletonVector(42, 2, 5)) == [0, 42, 0, 0, 0]
	@test size(SingletonVector(42, 5, 10_000)) == (10_000,)

	nan = SingletonVector("nan", 2, 5)
	@test nan[1] === 0 && nan[2] === "nan"
	@test eltype(nan) === Any
	@test collect(nan) == [0, "nan", 0, 0, 0]
end

@testset "@symbolicga" begin

	@basis 3
	x = (1, 0, 0)
	y = (0, 1, 0)
	z = (0, 0, 1)

	@test @symbolicga(3, (x=1, y=1), x + y, Tuple) == (1, 1, 0)
	@test @symbolicga(3, (x=1, y=1), x*y) == v12

	R = exp(π/4*v12)
	@assert grade(R) == 0:2:3
	@test @symbolicga(3, (x=1, R=0:2:3), grade(~R*x*R, 1)) ≈ v2

end