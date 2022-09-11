using Multivectors:
	EuclideanMetric,
	bits_of_grade

@testset "inverses" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)
	@test inv(v[1])v[1] == 1
	@test inv(5v[1]v[2]) == inv(v[2])inv(v[1])inv(5)

	@testset "dimension $dim" for dim in 0:5
		mixed_sig = (-1,+1,+1,+1,0)
		for sig in [
			EuclideanMetric(dim),
			mixed_sig[1:dim],
		]
			a = MixedMultivector{sig}(rand(2^dim))
			@test inv(a)*a ≈ 1
			@test a*inv(a) ≈ 1
		end
	end
end

@testset "exp" begin
	v = Blade{(1,1)}.(bits_of_grade(1, 2) .=> 1)
	@test exp(10000*2pi*v[1]v[2]) ≈ 1
	
	for dim in 1:5
		a = MixedMultivector{EuclideanMetric(dim)}(randn(2^dim))
		@test exp(a)exp(-a) ≈ 1
		@test inv(exp(a)) ≈ exp(-a) rtol=1e-4
	end
end