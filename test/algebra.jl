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

	o = Blade{(1,1)}.(1:3 .=> 0)
	@test o[1] == o[2] == o[3]

	v = Blade{(1,1,1)}(0b101 => 88.3)
	Ts = [Blade, Multivector, MixedMultivector]
	for T1 ∈ Ts, T2 ∈ Ts
		@test T1(v) == T2(v)
	end
end

@testset "≈" begin
	@test Blade{(1,1)}(0b10 => 1) ≈ Blade{(1,1)}(0b10 => 1 + eps())
	@test Blade{(1,1)}(0b10 => eps()) ≈ 0 atol=eps()

	# ≈ uses 2-norm on arrays, so atol ≥ √n*eps()
	@test Multivector{(1,1,1),1}(eps()rand(3)) ≈ 0 atol=√3*eps()
	@test Multivector{(1,1,1),0}([1]) ≈ 1 + eps() atol=eps()

	@test MixedMultivector{(1,1,1)}(1 .+ eps()rand(2^3)) ≈ MixedMultivector{(1,1,1)}(1 .+ eps()rand(2^3))

	v = Blade{(1,1,1)}.(0b101 .=> 10 .+ 1e-6rand(2))
	Ts = [Blade, Multivector, MixedMultivector]
	for T1 ∈ Ts, T2 ∈ Ts
		@test T1(v[1]) ≈ T2(v[2]) atol=1e-6
	end
end

@testset "scalar *" begin
	a = Blade{(1,1)}(0b01 => 10)

	for b in [a, Multivector(a), MixedMultivector(a)]
		@test 3a == a*3.0
		@test a/10 == 10\a
		@test -a == a/(-1.0)
	end

	@test a//5 === Blade{(1,1)}(0b01 => 2//1)
end

@testset "+" begin
	v = Blade{(0,0,0)}.(bits_of_grade(1, 3) .=> 1)
	bi = Blade{(0,0,0)}.(bits_of_grade(2, 3) .=> 1)

	@test v[1] + v[2] isa Multivector{(0,0,0),1}
	@test bi[1] + 2.5bi[2] isa Multivector{(0,0,0),2}
	@test v[1] + bi[2] isa MixedMultivector
	@test sum(Multivector.(v)) isa Multivector
	@test sum(MixedMultivector.(bi)) isa MixedMultivector
	@test v[1] - v[2] == -v[2] + v[1]

	@testset "with scalars" begin
		@test v[1] + 7.5 == 7.5 + v[1]
		@test (v[1] + v[2]) + 7 == v[1] + (v[2] + 7)
		@test (bi[1] + 8) - 8 == bi[1]
		@test 1 - v[1] == -(v[1] - 1)
		@test Multivector(v[1]) + 0.5 == v[1] + 1//2
	end
end



@testset "*" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)

	@test v[1]v[2] == -v[2]v[1]
	@test v[1]v[1] == -1
	@test (v[1] + 7v[3])v[2] == v[1]v[2] + 7v[3]v[2]
	@test (v[1]v[1] + v[2])v[2] == v[2]v[2] - v[2]
	@test v[1]v[2]v[3] == v[2]v[3]v[1] == v[3]v[1]v[2]
	@test 4 == 2*(v[2] + v[3])*(v[2] + v[3])

	Ts = [Blade, Multivector, MixedMultivector]
	for T1 in Ts, T2 in Ts
		@test v[1]v[2] == T1(v[1])T2(v[2]) == T2(v[1])T1(v[2])
	end
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

@testset "∧" begin
	v = Blade{(1,1,1)}.(bits_of_grade(1, 3) .=> 1)

	@test v[1]∧v[2] == -v[2]∧v[1] == v[1]v[2]
	@test (v[1] + v[2])∧v[2] == v[1]v[2]
	
end

@testset "^" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)
	
	@test v[1]^2 + 1 == 0
	@test v[1]^4003 == -v[1]

	@test (v[2] + v[3])^2 == 2
	@test (1 + v[1])^2 == 2v[1]

	a = 1 + 2v[2] + 3v[1]v[2]
	@test a^30 == ((a^2)^3)^5

	@test (v[1]v[2] + v[2]v[3])^10 == ((v[1]v[2])^2 + (v[2]v[3])^2)^5
end

@testset "reversion" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)
	
	@test reversion(v[1]) == v[1]
	@test reversion(v[1]v[2]) == v[2]v[1]
	@test reversion(3v[1]v[2] + v[3]v[4]) == 3v[2]v[1] + v[4]v[3]
	@test reversion(1 + 2v[1] + 3v[2]v[3] + 4v[1]v[2]v[3]) == 1 + 2v[1] + 3v[3]v[2] + 4v[3]v[2]v[1]

end