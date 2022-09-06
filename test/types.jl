using Multivectors:
	grade,
	signature

@testset "constructors and equality" begin
	@test zero(Blade{(1,1)}(0b11 => 42)) === Blade{(1,1),0,Int}(0b00, 0)
	@test one(Blade{(1,1)}(0b11 => 42.0)) === Blade{(1,1),0,Float64}(0b00, 1)

	@test zero(Multivector{(1,1,1),2,Vector{Int}}) == Multivector{(1,1,1),2}(zeros(Int, 3))
	@test zero(Multivector{(),0,Vector{Bool}}) == Multivector{(),0}([false])
end

@testset "Blade -> Multivector -> MixedMultivector" begin
	for a in [
		Blade{(1,1)}(0b00 => 42),
		Blade{(1,0,1)}(0b001 => 42.0),
		Blade{(1,1,1)}(0b111 => big(42)),
	]
		@test Multivector(a) isa Multivector{signature(a),grade(a),<:AbstractVector{eltype(a)}}
		@test MixedMultivector(a) isa MixedMultivector{signature(a),<:AbstractVector{eltype(a)}}
		@test MixedMultivector(Multivector(a)) isa MixedMultivector{signature(a),<:AbstractVector{eltype(a)}}
	end
end