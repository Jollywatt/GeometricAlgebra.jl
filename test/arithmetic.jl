using GeometricAlgebra

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
	@test a6/2 == 3a1
	@test a1 == 6.0\a6

	u = Multivector{(1,1,1),1}([1, 2, 4])
	v = Multivector{(1,1,1),1}([12, 24, 48])
	@test 12u == v
	@test 4(u*3) == v
	@test 120f0*u == BigInt(10)v
	@test u/2 == v/24
	@test 6\v == 2u

	A = MixedMultivector{(1,1)}([3, 0, 10, 7])
	B = MixedMultivector{(1,1)}([21, 0, 70, 49])
	@test 7A == B
	@test (7A)*2 == 2B
	@test (7 + 0im)A == (1//1)B
	@test B/7 == A
	@test 1\A == 7\B

end

@testset "addition, subtraction" begin
	sig = (1, 1, 1)
	x, y, z = basis(sig)
	@test x + y == Multivector{sig,1}([1, 1, 0])
	@test x + y == y + x
	@test (x + y) + z == y + x + z
	@test x + 2 == MixedMultivector{sig}([2, 1, 0, 0, 0, 0, 0, 0])
	@test x + 0 == x
	@test x + 3.0y == (3//1)y + x

	u = x + 2y
	@test -x == (-1)x
	@test -u == (-1)u
	@test +u == -(-u)
	@test x - y == x + (-y)
	@test u - 1 == -1 + u
	@test u - x == u - (0 + x)
end


@testset "geometric product" begin
	t, x, y, z = basis((-1, 1, 1, 1))
	@testset "Blade" begin
		@test x*x == y*y == z*z == 1
		@test t*t == -1
		@test x*y == -y*x
		@test (x*y)*(x*y) == -1
		@test (y*z)*(t*x) == t*x*y*z
		@test x*x*x == x
	end

	@testset "Multivector" begin
		u = 3x + 4y
		@test u*x == 3 + 4y*x
		@test u*u == 25
		@test (t + x)*(t + x) == 0
	end

	@testset "MixedMultivector" begin
		@test (1 + x)*(1 + y) == 1 + x + y + x*y
		@test (1 + x)*(1 - x) == 0
		@test (2 + x*y)*(y + x*y) == 2y + 2x*y + x - 1
	end
end


@testset "multiplicative inverses" begin
	t, x, y, z = basis((-1, 1, 1, 1))
    @test inv(x) == x
    @test inv(t) == -t
    @test inv(x*y) == y*x
    @test x/x == y\y == 1
    @test (t*x*y*z)/(t*z) == t*x*y*z*t*z == x*y

    @test inv(3x + 4y) == (3x + 4y)/25
    cases = [
    	x + 2y + 3z,
    	t*x + y*z, # simplest multivector with non-scalar square
    	1 + 2x,
    	3 + x*y*z,
    	1 + t + x*y
    ]
    @testset "inv($a)" for a ∈ cases
    	@test inv(a)*a ≈ 1
    	@test a*inv(a) ≈ 1
    end
end


@testset "integer powers" begin
	t, x, y, z = basis((-1, 1, 1, 1))
	
	@testset "power_with_scalar_square()" begin
		using GeometricAlgebra: power_with_scalar_square
		pwss = power_with_scalar_square

		@testset "a = $a" for a ∈ float.([x, 2x*y, x*t])
			@testset "a^p, p = $p" for p ∈ [-4:4; 13; 14; 171; -333; 10042]
				a² = a*a
				@assert isscalar(a²)
				@test pwss(a, scalar(a²), p) == prod(fill(p > 0 ? a : inv(a), abs(p)))
			end
		end
	    
	end
end



@testset "duality operations" begin
	t, x, y, z = basis((-1, 1, 1, 1))
	@test reversion(x*y) == y*x
	@test reversion(42) == 42
	@test ~x == x
	@test ~(x*y + y*z) == y*x + z*y
	@test ~(1 + t*x*y) == 1 - t*x*y
	for a ∈ [x*z + y*z, 1 + x + x*y + x*y*z + t*x*y*z]
		@test ~(~a) == a
	end

	@test involute(4) == 4
	@test involute(4x) == -4x
	@test involute(4x*y) == 4x*y
	@test involute(1 + 4x*y + t) == 1 + 4x*y - t
end


@testset "grade operations" begin
	t, x, y, z = basis((-1, 1, 1, 1))
    @test grade(42) == 0
    @test grade(x) == 1
    @test grade(x*y) == 2
    @test grade(x + y) == 1

    @test grade(42, 1) == 0
    @test grade(x, 1) == x
    @test grade(x, 10) == 0
    @test grade(x + y, 1) == x + y
    @test grade(x + y, 0) == 0
    @test grade(1 + x + y*z, 2) == y*z
end