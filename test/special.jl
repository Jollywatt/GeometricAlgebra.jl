using Multivectors:
	bits_of_grade

@testset "inverses" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1,4) .=> 1)

	@testset "scalar a²" begin
		for a in [
			v[2],
			5v[1]v[2],
			v[2]v[3]v[4]//10,
			v[2] + 4v[4],
			2.0v[1] + π*v[2],
			v[2] + 2v[2]v[3]
		]
			@test a*inv(a) == 1 == inv(a)*a
			@test a^3/a == a^2
		end
	end

	@testset "scalar aã" begin
		for a in [
			v[1]v[2] + 3
		]
			@test a*inv(a) == 1 == inv(a)*a
		end
	end
end