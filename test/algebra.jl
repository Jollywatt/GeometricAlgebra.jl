@testset "==" begin
	for sig in [(1,1), (-1,0,+1)], T in [Bool, Int, Float64], k in 1:3
		@test zero(Blade{sig,k,T}) == Blade{sig}(0 => zero(T))
		@test one(Blade{sig,k,T}) == Blade{sig}(0 => one(T))

		@test zero(Multivector{sig,k,Vector{T}}) == Multivector{sig,k}(zeros(T, binomial(dimension(sig), k)))
		@test zero(MixedMultivector{sig,Vector{T}}) == MixedMultivector{sig}(zeros(T, 2^dimension(sig)))
	end

	@test Blade{(1,1)}(0b00 => 0) == Blade{(1,1)}(0b11 => 0)
	@test zero(Multivector{(1,1,1),0,Vector{Int}}) == zero(Multivector{(1,1,1),3,Vector{Int}})
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