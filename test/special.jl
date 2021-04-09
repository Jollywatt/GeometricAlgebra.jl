using GeometricAlgebra

@testset "trig" begin
	t, x, y, z = basis((-1, 1, 1, 1))

	@testset "scalar square" begin
		@test exp(+x) ≈ cosh(1) + sinh(1)x
		@test exp(0y) ≈ 1
		@test exp(-z) ≈ cosh(1) - sinh(1)z

		@test exp(x*y) ≈ cos(1) + sin(1)x*y

		for a ∈ (x, 2y*x, x + y, x*y + 3y*z)
			@test tan(a) ≈ sin(a)/cos(a)
			@test tanh(a) ≈ sinh(a)/cosh(a)
			@test cot(a) ≈ cos(a)/sin(a)
			@test coth(a) ≈ cosh(a)/sinh(a)

			@test cos(a)^2 + sin(a)^2 ≈ 1
			@test 1 + tan(a)^2 ≈ sec(a)^2
			@test cot(a)^2 + 1 ≈ csc(a)^2
			
			@test cosh(a)^2 - sinh(a)^2 ≈ 1
			@test 1 - tanh(a)^2 ≈ sech(a)^2
			@test coth(a)^2 - 1 ≈ csch(a)^2
		end
	end
end