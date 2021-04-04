using GeometricAlgebra

@testset "constructors" begin

	@test iszero(zero(Blade{(1,1,1),2,0b101,Float64}))
	@test iszero(zero(Blade{(1,1,1),0,0b000,BigInt}))
	@test iszero(zero(Multivector{(1,1,1),1,Vector{Int}}))
	@test iszero(zero(MixedMultivector{(1,1,1),Vector{Int}}))

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
end

@testset "scalar multiplication" begin
	a1 = Blade{(1,1)}(1, 0b01)
	a6 = Blade{(1,1)}(6, 0b01)
	@test 0a1 == 0a6
	@test a1*6 == a6
	@test 6*a1 == a6
	@test (3*a1)*2 == a6
	@test a6*2 == 3(2a1)*2

	@test 3.0a1 == 0.5a6
	@test 6im*a1 == im*a6

	u = Multivector{(1,1,1),1}([1, 2, 4])
	v = Multivector{(1,1,1),1}([12, 24, 48])
	@test 12u == v
	@test 4(u*3) == v
	@test 120f0*u == BigInt(10)v

	A = MixedMultivector{(1,1)}([3, 0, 10, 7])
	B = MixedMultivector{(1,1)}([21, 0, 70, 49])
	@test 7A == B
	@test (7A)*2 == 2B
	@test (7 + 0im)A == (1//1)B

end