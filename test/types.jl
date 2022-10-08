using GeometricAlgebra:
	grade,
	bits_of_grade,
	isscalar,
	signature,
	largest_type,
	MetricWithStorage
using GeometricAlgebra.StaticArrays
using GeometricAlgebra.SparseArrays

@testset "constructors" begin
	@test zero(Blade{2}(0b11 => 42)) === Blade{2,0,Int}(0b00, 0)
	@test one(Blade{(1,1)}(0b11 => 42.0)) === Blade{(1,1),0,Float64}(0b00, 1)

	for T in [
		Int,
		Blade{(1,1),0,Bool},
		KVector{(1,1,1),2,Vector{Int}},
		KVector{(),0,Vector{Bool}},
		KVector{(0,1),2,Vector{Float64}},
		Multivector{(-1,+1,+1,+1),Vector{Float64}},
	]
		@test zero(T) isa T
		@test zero(zero(T)) isa T
		@test iszero(zero(T))
		@test isone(one(T))
		@test isscalar(zero(T))
	end
end

@testset "Blade -> KVector -> Multivector" begin
	for a in [
		Blade{(1,1)}(0b00 => 42),
		Blade{(1,0,1)}(0b001 => 42.0),
		Blade{3}(0b111 => big(42)),
	]
		@test KVector(a) isa KVector{signature(a),grade(a),<:AbstractVector{eltype(a)}}
		@test Multivector(a) isa Multivector{signature(a),<:AbstractVector{eltype(a)}}
		@test Multivector(KVector(a)) isa Multivector{signature(a),<:AbstractVector{eltype(a)}}
	end

	x = Blade{(1,1)}(0b01 => 5)

	Ts = [Blade, KVector, Multivector]
	for T1 in Ts, T2 in Ts
		a, b = T1(x), T2(x)
		@test largest_type(a, b) == largest_type(b, a)
		@test largest_type(a, b)(a) isa largest_type(a, b)
	end

	for (M1, M2) in [
		Blade => KVector,
		KVector => KVector,
		Blade => Multivector,
		KVector => Multivector,
		Multivector => Multivector,
	], T in [Int, Float64, ComplexF32]
		@test eltype(M2(M1(x), T)) === T
	end

	@test_throws ErrorException length(x)
	@test ncomponents(KVector(x)) == 2
	@test ncomponents(Multivector(x)) == 4
end

@testset "eltype" begin
	for T in [Int, Float64, Complex{Float16}, Rational{Int128}, Bool]
		a = Blade{(1,1)}(0b11 => one(T))

		@test eltype(a) === T
		@test eltype(KVector(a)) === T
		@test eltype(Multivector(a)) === T
	end
end

@testset "different storage types" begin
	@testset "$S" for S in [SVector, MVector, SparseVector]
		sig = MetricWithStorage{(1,1,1),S}()
		v = basis(sig)

		@test KVector(v[1]) isa KVector{sig,1,<:S}
		@test Multivector(v[1]) isa Multivector{sig,<:S}

		for T in [Blade, KVector, Multivector]
			@test isone(one(T(v[1])))
		end

		for (M1, M2) in [
			Blade => KVector,
			KVector => KVector,
			Blade => Multivector,
			KVector => Multivector,
			Multivector => Multivector,
		], T in [Int, Float64, ComplexF32]
			@test eltype(M2(M1(v[1]), T)) === T
		end
	end
end