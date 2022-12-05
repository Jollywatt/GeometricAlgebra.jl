using GeometricAlgebra:
	zeroslike,
	oneslike,
	symbolic_multivector
using GeometricAlgebra.SymbolicUtils


@testset "symbolic components" begin
	@test zeroslike(Vector{Int}, 1) == [0]
	@test zeroslike(Vector{Any}, 2) == [0, 0]
	@test zeroslike(Vector{SymbolicUtils.Symbolic}, 3) == [0, 0, 0]
	@test zeroslike(Vector{Union{Int,SymbolicUtils.Symbolic}}, 4) == [0, 0, 0, 0]


	b = BasisBlade{3}(0b101 => SymbolicUtils.Sym{Real}(:b))
	kv = symbolic_multivector(Multivector{3,1,Vector{Int}}, :kv)
	mv = symbolic_multivector(Multivector{3,0:3,Vector{Int}}, :mv)

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