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


	b = Blade{(1,1,1)}(0b101 => SymbolicUtils.Sym{Real}(:b))
	m = symbolic_multivector(Multivector{(1,1,1),1,Vector{Int}}, :m)
	mm = symbolic_multivector(MixedMultivector{(1,1,1),Vector{Int}}, :mm)

	for a in [b, m, mm]
		@test !iszero(a)
		@test !isone(a)

		@test iszero(zero(a))
		@test isone(one(a))

		@test MixedMultivector(a) isa MixedMultivector
	end

	@testset "symbolic algebra" begin
		@test b + 2m + π*mm isa MixedMultivector
		@test grade(b^2) == 0
		@test m*5mm isa MixedMultivector
		@test 1 + m isa MixedMultivector

		@test m.^(0:3) isa Vector{<:AbstractMultivector}

		@test iszero(m∧m)
		@test b∧m∧mm isa MixedMultivector
	end

end