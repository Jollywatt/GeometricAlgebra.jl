using GeometricAlgebra:
	bits_of_grade,
	scalar

@testset "inverses" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)
	@test inv(v[1])v[1] == 1
	@test inv(5v[1]v[2]) == inv(v[2])inv(v[1])inv(5)

	u = sum(v)
	@test 1/u == u\1 == u/scalar(u^2) == scalar(u^2)\u

	@testset "dimension $dim" for dim in 0:5
		mixed_sig = (-1,+1,+1,+1,0)
		for sig in [
			dim,
			mixed_sig[1:dim],
		], trials in 1:5
			a = Multivector{sig}(rand(2^dim))
			@test inv(a)*a ≈ 1 rtol=1e-6
			@test 1 ≈ inv(a)*a rtol=1e-6
		end
	end
end

@testset "exp" begin
	v = Blade{(1,1)}.(bits_of_grade(1, 2) .=> 1)
	@test exp(10000*2pi*v[1]v[2]) ≈ 1

	@test exp(0v[1]) == 1
	@test exp(3v[1]) ≈ cosh(3) + sinh(3)v[1]
	@test exp(3v[1]v[2]) ≈ cos(3) + sin(3)v[1]v[2]
	
	for dim in 1:5, trials in 1:5
		a = Multivector{dim}(big.(randn(2^dim)))
		@test exp(a)exp(-a) ≈ 1
		@test inv(exp(a)) ≈ exp(-a)
	end
end

@testset "matrix method fallbacks" begin
	v = basis(3)

	a = 1 + v[1] + 2v[1]v[3]
	@test log(exp(a)) ≈ a

	@test exp(log(a)/2) ≈ sqrt(a)
end

@testset "trig identities" begin
	for sig in [3, (0, -1, 1), 4]
		v = basis(sig)
		for a in [5v[1]v[2], v[1]v[2] + v[2]v[3], rand(length(v))'v]
			@test exp(a) ≈ cosh(a) + sinh(a)
			@test tan(a) ≈ sin(a)/cos(a)
			@test sin(a) ≈ sin(asin(sin(a)))
		end
	end
end