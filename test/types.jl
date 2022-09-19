using GeometricAlgebra:
	grade,
	isscalar,
	signature,
	largest_type,
	with_eltype,
	MetricWithStorage
using GeometricAlgebra.StaticArrays	

@testset "constructors" begin
	@test zero(Blade{(1,1)}(0b11 => 42)) === Blade{(1,1),0,Int}(0b00, 0)
	@test one(Blade{(1,1)}(0b11 => 42.0)) === Blade{(1,1),0,Float64}(0b00, 1)

	for T in [
		Int,
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
		@test isscalar(zero(T))
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

	for (M1, M2) in [
		Blade => Multivector,
		Multivector => Multivector,
		Blade => MixedMultivector,
		Multivector => MixedMultivector,
		MixedMultivector => MixedMultivector,
	], T in [Int, Float64, ComplexF32]
		@test eltype(M2(M1(x), T)) === T
	end

	@test_throws ErrorException length(x)
	@test ncomponents(Multivector(x)) == 2
	@test ncomponents(MixedMultivector(x)) == 4
end

@testset "eltype" begin
	for T in [Int, Float64, Complex{Float16}, Rational{Int128}, Bool]
		a = Blade{(1,1)}(0b11 => one(T))

		@test eltype(a) === T
		@test eltype(Multivector(a)) === T
		@test eltype(MixedMultivector(a)) === T
	end
end

@testset "different storage types" begin
	for S in [SVector, MVector]
		sig = MetricWithStorage{(1,1,1),S}()
		v = Blade{sig}.(bits_of_grade(1, 3) .=> 1)

		@test Multivector(v[1]) isa Multivector{sig,1,<:S}
		@test MixedMultivector(v[1]) isa MixedMultivector{sig,<:S}

		for T in [Blade, Multivector, MixedMultivector]
			@test isone(one(T(v[1])))
		end

		for (M1, M2) in [
			Blade => Multivector,
			Multivector => Multivector,
			Blade => MixedMultivector,
			Multivector => MixedMultivector,
			MixedMultivector => MixedMultivector,
		], T in [Int, Float64, ComplexF32]
			@test eltype(M2(M1(v[1]), T)) === T
		end
	end
end