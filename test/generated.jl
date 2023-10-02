using GeometricAlgebra:
	zeroslike,
	make_symbolic
import GeometricAlgebra.SymbolicUtils
using GeometricAlgebra.StaticArrays:
	MVector,
	SVector


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

@testset "@symbolicga" begin

	@basis 3
	x = (1, 0, 0)
	y = (0, 1, 0)
	z = (0, 0, 1)

	@test @symbolicga(3, (x=1, y=1), x + y, SVector) === SVector(1, 1, 0)
	@test @symbolicga(3, (x=1, y=1), x*y) == v12

	R = exp(π/4*v12)
	@assert grade(R) == 0:2:3
	@test @symbolicga(3, (x=1, R=0:2:3), grade(~R*x*R, 1)) ≈ v2

	joinpoints(a::Tuple, b::Tuple) = @symbolicga 3 (a=1, b=1) a∧b Tuple
	meetlines(a::Tuple, b::Tuple) = @symbolicga 3 (a=2, b=2) a∨b Tuple
	from_homogeneous((x, y, z)) = (x/z, y/z)

	l1 = @inferred joinpoints((1, 0, 1), (0, 1, 1)) # line y = 1 - x 
	l2 = @inferred joinpoints((0, 0, 1), (1, 1, 1)) # line y = x
	@test from_homogeneous(meetlines(l1, l2)) == (0.5, 0.5)

end