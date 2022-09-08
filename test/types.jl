using Multivectors:
	grade,
	signature,
	largest_type

@testset "constructors" begin
	@test zero(Blade{(1,1)}(0b11 => 42)) === Blade{(1,1),0,Int}(0b00, 0)
	@test one(Blade{(1,1)}(0b11 => 42.0)) === Blade{(1,1),0,Float64}(0b00, 1)

	for T in [
		Blade{(1,1),0,Bool},
		Multivector{(1,1,1),2,Vector{Int}},
		Multivector{(),0,Vector{Bool}},
		Multivector{(0,1),2,Vector{Float64}},
		MixedMultivector{(-1,+1,+1,+1),Vector{Float64}},
	]
		@test zero(T) isa T
		@test zero(zero(T)) isa T
		@test iszero(zero(T))

		@test isone(one(T))
	end
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

	x = Blade{(1,1)}(0b01 => 5)
	Ts = [Blade, Multivector, MixedMultivector]
	for T1 in Ts, T2 in Ts
		a, b = T1(x), T2(x)
		@test largest_type(a, b) == largest_type(b, a)
		@test largest_type(a, b)(a) isa largest_type(a, b)
	end
end