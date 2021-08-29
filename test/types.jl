using GeometricAlgebra
using GeometricAlgebra: StaticArrays.SVector

@testset "constructors" begin
	# @test iszero(zero(Blade{(1,1,1),2,0b101,Float64}))
	@test iszero(zero(Blade{(1,1,1)}(big(42), 0b000)))
	@test iszero(zero(Multivector{(1,1,1),1,Vector{Int}}))
	@test iszero(zero(MixedMultivector{(1,1,1),Vector{Int}}))

	@test first(basis((1,1,1))) == Blade{(1,1,1)}(1, 0b001)
	@test  last(basis((1,1,1))) == Blade{(1,1,1)}(1, 0b100)
end


@testset "equality" begin
	a = Blade{(1,1)}(42, 0b01)
	@test a == Blade{(1,1)}(42, 0b01)
	@test a == Blade{(1,1)}(42.0, 0b01)
	@test a == Blade{(1,1)}(BigFloat(42), 0b01)
	@test a != Blade{(1,1)}(-1, 0b01)
	@test a != Blade{(1,1)}(42, 0b10)
	@test a != Blade{(1,1)}(42, 0b11)
	@test a != Blade{(1,1,1)}(42, 0b01)

	a = Multivector{(1,1,1),1}([3, 5, 7])
	@test a == Multivector{(1,1,1),1}(Float16[3, 5, 7])
	@test a == Multivector{(1,1,1),1}(Number[3, 5, 7 + 0im])
	@test a != Multivector{(1,1,1),1}([3, -1, 7])
	@test a != Multivector{(1,1,1),2}([3, 5, 7])
	@test a != Multivector{(-1,1,1),1}([3, 5, 7])

	a = MixedMultivector{(1,1)}([1, 10, 20, 100])
	@test a == MixedMultivector{(1,1)}(BigInt[1, 10, 20, 100])
	@test a != MixedMultivector{(1,1)}([1, 10, -1, 100])
	@test a != MixedMultivector{(-1,1)}([1, 10, 20, 100])

	# equal if both zero even if grade differs
	@test Blade{(1,1)}(0, 0b01) == Blade{(1,1)}(0, 0b10)
	@test Multivector{(1,1),1}([0, 0]) == Multivector{(1,1),0}([0])
	@test Blade{(1,1)}(0, 0b01) != Blade{(1,1)}(1, 0b10)

	@test Blade{(1,1)}(0, 0b01) == 0
	@test Blade{(1,1)}(1, 0b00) == 1
	@test Blade{(1,1)}(2, 0b00) != 1
	@test Blade{(1,1)}(1, 0b01) != 1
end


# TODO: test getcomponent and Base.getindex


@testset "best_type()" begin
	x = Blade{(1,1),1,Int}(42, 0b01)

	@testset "given $label" for (label, a) ∈ ["intances" => x, "types" => typeof(x)]

		@test best_type(Blade, a) === Blade{(1,1),k,Int} where k
		@test best_type(Multivector, a) === Multivector{(1,1),k,Vector{Int}} where k
		@test best_type(MixedMultivector, a) === MixedMultivector{(1,1),Vector{Int}}

		@test best_type(Blade, a; grade=Val(2)) === Blade{(1,1),2,Int}
		@test best_type(Blade, a; promote_eltype_with=BigInt) === Blade{(1,1),k,BigInt} where k
		@test best_type(Blade, a; grade=Val(0), promote_eltype_with=Bool) === Blade{(1,1),0,Int}
		@test best_type(Multivector, a; grade=Val(2)) === Multivector{(1,1),2,Vector{Int}}
		@test best_type(Multivector, a; promote_eltype_with=Complex{Int}) === Multivector{(1,1),k,Vector{Complex{Int}}} where k
		@test best_type(Multivector, a; grade=Val(1), promote_eltype_with=Float32) === Multivector{(1,1),1,Vector{Float32}}
		@test best_type(MixedMultivector, a; promote_eltype_with=Float32) === MixedMultivector{(1,1),Vector{Float32}}

		# single-argument versions preserve parameters
		@test best_type(x) == (x isa Type ? x : typeof(x))
		@test best_type(a; promote_eltype_with=BigFloat, grade=Val(0)) === Blade{(1,1),0,BigFloat}
	end

	y = Multivector{(1,1),1,Vector{Float64}}(0:1)
	# @test best_type(y) == typeof(y)
	# @test best_type(Multivector, x, y) == Multivector{(1,1),k,Vector{Float64}} where k

	# @test_throws Exception best_type(MixedMultivector, MixedMultivector{(1,1,1),Vector{Int}}, Multivector{(1,1),0,Vector{Int}})
end


@testset "conversion" begin
	scalar_types = [
		Int,
		Int128,
		BigInt,
		Float64,
		BigFloat,
		Complex{Int},
		Complex{Float64},
		Rational{Int8},
	]
	@testset "eltype conversions" begin
		blades = Blade{(1,)}.(one.(scalar_types), 0b1)
		@testset "$(eltype(a)) -> $(eltype(b))" for a ∈ blades, b ∈ blades
			@test convert(typeof(b), a) == b
		end
	end

	x, y = basis((1, 1))
	@test float(x) isa Blade{(1, 1),k,<:AbstractFloat} where k
	@test big(x) isa Blade{(1, 1),k,BigInt} where k
	@test complex(x) isa Blade{(1, 1),k,Complex{Int}} where k

	fns = float, big, complex, real
	objs = x*y, x + y, x*y - 1
	@testset "$(GeometricAlgebra.multivectortype(a))" for a ∈ objs
		@testset "$fn()" for fn ∈ fns
			@test typeof(fn(a)) == fn(typeof(a))
			@test typeof(fn(complex(a))) == fn(typeof(complex(a)))
		end
		@test real(a + 2im*a) == a
		@test imag(2a + im*a) == a
	end

	@test Multivector(x)::Multivector == x
	@test MixedMultivector(x)::MixedMultivector == x
	@test MixedMultivector(x + 1) == x + 1
end


@testset "promotion" begin
	x, y, z = basis((1, 1, 1))
	@test ==(promote(2x, x + x)...)
	@test ==(promote(x, x + 0)...)
	@test ==(promote(x + y, x + y + 0)...)

	@test ==(promote(1, x^2)...)
	@test ==(promote(2, (x + y)^2)...)
	@test ==(promote(1, (x + 0)^2)...)

	# objects are promoted when stored in vectors
	@test first([1, x]) isa Blade
	@test first([2, x + y]) isa Multivector
	@test first([3, z + 1]) isa MixedMultivector
	@test sum([1, x, x + y, z + 1]) == 2 + 2x + y + z

	@test prod([x, x*y, y*x*z]) == x*z
end


@testset "isapprox()" begin
	x, y, z = basis((1, 1, 1))
	@test x^2 ≈ 1
	@test sin(1f0)x ≈ sin(1)x
	@test (x + y)^2 ≈ 2
	@test (1 + x*y)^4 ≈ -4

	@test ((x + y)/√2)^Int(1e10) ≈ 1 rtol=1e-4

	# TODO: test zero-equalities
end

@testset "storagetype" begin
	sig = (1, 1, 1)
	b = basis(sig)
	B = convert.(Multivector{sig,1,SVector{k,Int} where k}, b)

	@test b[1]b[2] isa Blade{sig,2}
	@test B[1]B[2] isa MixedMultivector{sig,<:SVector}

	
	# we want operations with blades to preserve the storage type 
	isstatic(a) = GeometricAlgebra.storagetype(a) <: SVector
	@test isstatic(B[1] + b[2])
	@test isstatic(B[1] - b[2])
	@test isstatic(B[1]b[2])
	@test isstatic(B[1]/b[2])
	@test isstatic((B[1] + 1)∧b[2])
	@test isstatic((B[1] + 1)*b[2])
end