using Multivectors:
	bits_of_grade

@testset "==" begin
	for sig in [(1,1), (-1,0,+1)], T in [Bool, Int, Float64], k in 1:3
		@test zero(Blade{sig,k,T}) == Blade{sig}(0 => zero(T))
		@test one(Blade{sig,k,T}) == Blade{sig}(0 => one(T))

		@test zero(Multivector{sig,k,Vector{T}}) == Multivector{sig,k}(zeros(T, binomial(dimension(sig), k)))
		@test zero(MixedMultivector{sig,Vector{T}}) == MixedMultivector{sig}(zeros(T, 2^dimension(sig)))
	end

	@test Blade{(1,1)}(0b00 => 0) == Blade{(1,1)}(0b11 => 0)
	@test zero(Multivector{(1,1,1),0,Vector{Int}}) == zero(Multivector{(1,1,1),3,Vector{Int}})

	v = Blade{(1,1)}.(1:3 .=> 0)
	@test v[1] == v[2] == v[3]
end

@testset "scalar *" begin
	a = Blade{(1,1)}(0b01 => 10)

	for b in [a, Multivector(a), MixedMultivector(a)]
		@test 3a == a*3
		@test a/10 == 10\a
		@test -a == a/(-1)
	end

	@test a//5 === Blade{(1,1)}(0b01 => 2//1)
end

@testset "+" begin
	v = Blade{(0,0,0)}.(bits_of_grade(1, 3) .=> 1)
	bi = Blade{(0,0,0)}.(bits_of_grade(2, 3) .=> 1)

	@test v[1] + v[2] isa Multivector{(0,0,0),1}
	@test bi[1] + bi[2] isa Multivector{(0,0,0),2}
	@test v[1] + bi[2] isa MixedMultivector
	@test sum(Multivector.(v)) isa Multivector
	@test sum(MixedMultivector.(bi)) isa MixedMultivector
	@test v[1] - v[2] == -v[2] + v[1]

	@testset "with scalars" begin
		@test v[1] + 7 == 7 + v[1]
		@test (v[1] + v[2]) + 7 == v[1] + (v[2] + 7)
		@test (bi[1] + 8) - 8 == bi[1]
		@test 1 - v[1] == -(v[1] - 1)
	end
end



@testset "*" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)

	@test v[1]v[2] == -v[2]v[1]
	@test v[1]v[1] == -1
	@test (v[1] + 7v[3])v[2] == v[1]v[2] + 7v[3]v[2]
	@test (v[1]v[1] + v[2])v[2] == v[2]v[2] - v[2]
	@test v[1]v[2]v[3] == v[2]v[3]v[1] == v[3]v[1]v[2]

	Ts = [Blade, Multivector, MixedMultivector]
	for T in Ts
		@test v[1]v[2] == T(v[1])v[2] == v[1]T(v[2]) == T(v[1])T(v[2])
	end

	@test Multivector(v[1])MixedMultivector(v[2]v[3]) == v[1]v[2]v[3]
end

@testset "scalar product" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)

	@test scalar_prod(v[1], v[1]) == -1
	@test scalar_prod(v[1], v[2]) == 0
	@test typeof(scalar_prod(v[1], v[1])) === typeof(scalar_prod(v[1], v[2]))

	@test scalar_prod(v[2] + v[3], v[2] + v[3]) == 2
	@test scalar_prod(v[1] + v[2], v[3]) == 0
	@test scalar_prod(10 + 2v[3], 20 + v[3]) == 202
end

@testset "^" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)
	
	@test v[1]^2 == -1
	@test v[1]^4003 == -v[1]

	@test (v[2] + v[3])^2 == 2
	@test (1 + v[1])^2 == 2v[1]

	a = 1 + 2v[2] + 3v[1]v[2]
	@test a^30 == ((a^2)^3)^5
end