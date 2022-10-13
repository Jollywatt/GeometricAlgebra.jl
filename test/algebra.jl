using GeometricAlgebra:
	bits_of_grade,
	MetricWithStorage,
	isscalar,
	scalar,
	unit_pseudoscalar,
	isapproxzero

using GeometricAlgebra.StaticArrays
using GeometricAlgebra.SparseArrays

@testset "==" begin
	for sig in [(1,1), (-1,0,+1)], T in [Bool, Int, Float64], k in 1:3
		@test zero(Blade{sig,k,T}) == Blade{sig}(0 => zero(T))
		@test one(Blade{sig,k,T}) == Blade{sig}(0 => one(T))

		@test zero(KVector{sig,k,Vector{T}}) == KVector{sig,k}(zeros(T, binomial(dimension(sig), k)))
		@test zero(Multivector{sig,Vector{T}}) == Multivector{sig}(zeros(T, 2^dimension(sig)))
	end

	@test Blade{(1,1)}(0b00 => 0) == Blade{(1,1)}(0b11 => 0)
	@test zero(KVector{(1,1,1),0,Vector{Int}}) == zero(KVector{(1,1,1),3,Vector{Int}})

	o = Blade{(1,1)}.(1:3 .=> 0)
	@test o[1] == o[2] == o[3]

	v = Blade{(1,1,1)}(0b101 => 88.3)
	Ts = [Blade, KVector, Multivector]
	for T1 ∈ Ts, T2 ∈ Ts
		@test T1(v) == T2(v)
	end

end

@testset "≈" begin
	@test Blade{(1,1)}(0b10 => 1) ≈ Blade{(1,1)}(0b10 => 1 + eps())
	@test Blade{(1,1)}(0b10 => eps()) ≈ 0 atol=eps()

	# ≈ uses 2-norm on arrays, so atol ≥ √n*eps()
	@test KVector{(1,1,1),1}(eps()rand(3)) ≈ 0 atol=√3*eps()
	@test KVector{(1,1,1),0}([1]) ≈ 1 + eps() atol=eps()

	@test Multivector{(1,1,1)}(1 .+ eps()rand(2^3)) ≈ Multivector{(1,1,1)}(1 .+ eps()rand(2^3))

	v = Blade{(1,1,1)}.(0b101 .=> 10 .+ 1e-6rand(2))
	Ts = [Blade, KVector, Multivector]
	for T1 ∈ Ts, T2 ∈ Ts
		@test T1(v[1]) ≈ T2(v[2]) atol=1e-6
	end

	v = Blade{(1,1,1)}(0b110 => 1)
	@test isapproxzero(eps()v; atol=eps())
	@test isapproxzero(eps()v + 0; atol=eps())
end

@testset "grade projections" begin
	@basis 3

	@test isscalar(scalar(v1))
	@test scalar(v1 + v2) == 0
	@test scalar(v1 + 42) == 42

	@test grade(v1 + v12 + v2, 1) == v1 + v2
	@test grade(v1 + v2, 0) == 0

	u = v1 + 2v2 + 3v3
	@test grade(u, grade(u)) === u
end

@testset "scalar *" begin
	a = Blade{(1,1)}(0b01 => 10)

	for b in [a, KVector(a), Multivector(a)]
		@test 3a == a*3.0
		@test a/10 == 10\a
		@test -a == a/(-1.0)
	end

	@test a//5 === Blade{(1,1)}(0b01 => 2//1)
end

@testset "+" begin
	v = Blade{(0,0,0)}.(bits_of_grade(1, 3) .=> 1)
	bi = Blade{(0,0,0)}.(bits_of_grade(2, 3) .=> 1)

	@test v[1] + v[2] isa KVector{(0,0,0),1}
	@test bi[1] + 2.5bi[2] isa KVector{(0,0,0),2}
	@test v[1] + bi[2] isa Multivector
	@test sum(KVector.(v)) isa KVector
	@test sum(Multivector.(bi)) isa Multivector
	@test v[1] - v[2] == -v[2] + v[1]

	@testset "with scalars" begin
		@test v[1] + 7.5 == 7.5 + v[1]
		@test (v[1] + v[2]) + 7 == v[1] + (v[2] + 7)
		@test (bi[1] + 8) - 8 == bi[1]
		@test 1 - v[1] == -(v[1] - 1)
		@test KVector(v[1]) + 0.5 == v[1] + 1//2
	end

	@test (1:3)'v == 1v[1] + 2v[2] + 3v[3]
end

@testset "*" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)

	@test v[1]v[2] == -v[2]v[1]
	@test v[1]v[1] == -1
	@test (v[1] + 7v[3])v[2] == v[1]v[2] + 7v[3]v[2]
	@test (v[1]v[1] + v[2])v[2] == v[2]v[2] - v[2]
	@test v[1]v[2]v[3] == v[2]v[3]v[1] == v[3]v[1]v[2]
	@test 4 == 2*(v[2] + v[3])*(v[2] + v[3])

	@test GeometricAlgebra.geometric_prod(v[1], 2) == 2v[1]

	Ts = [Blade, KVector, Multivector]
	for T1 in Ts, T2 in Ts
		@test v[1]v[2] == T1(v[1])T2(v[2]) == T2(v[1])T1(v[2])
	end
end

@testset "scalar product" begin
	v = Blade{(-1,+1,+1,+1)}.(bits_of_grade(1, 4) .=> 1)

	@test v[1] ⊙ v[1] == -1
	@test v[1] ⊙ v[2] == 0
	@test v[1] ⊙ (1 + v[1]v[2]) isa Real

	@test scalar_prod(v[2] + v[3], v[2] + v[3]) == 2
	@test scalar_prod(v[1] + v[2], v[3]) == 0
	@test scalar_prod(10 + 2v[3], 20 + v[3]) == 202
end

@testset "∧" begin
	v = Blade{(1,1,1)}.(bits_of_grade(1, 3) .=> 1)

	@test v[1]∧v[2] == -v[2]∧v[1] == v[1]v[2]
	@test (v[1] + v[2])∧v[2] == v[1]v[2]

	@test v[1]∧10 == 2∧v[1]∧5 == 10v[1]
end

@testset "⋅" begin
	v = basis(3)
	@test v[1]⋅v[1] == v[2]⋅v[2] == 1
	@test v[1]⋅v[2] == 0
	@test (v[1]v[2])⋅v[2] == v[2]⋅(v[2]v[1]) == v[1]
	@test (1 + v[1] + 2v[2]v[3])⋅v[3] == v[3] + 2v[2]
	@test v[1]⋅8 == 2⋅v[1]⋅(2⋅2)

	v = basis("-+++")
	@test (1 + v[1])⋅(10 + 5v[1]) == 10 + (5 + 10)v[1] - 5
end

@testset "⨼, ⨽" begin
	@basisall 3

	@test v123⨽v32 == v1 == ~(v23⨼v321)
	@test iszero(7⨽(v1 + v2))
	@test v2⨽(10 + v123) == 10v2
	@test 2⨼2 == 4
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

@testset "different storage types" begin
	for S in [SVector, MVector, SparseVector]
		sig = MetricWithStorage{(1,1,1),S}()
		v = basis(sig)

		@test v[1]v[2] == ~(v[2]v[1])
		@test 5.5Multivector(v[1]) + v[3] == v[3] + 5.5v[1]
		@test 1 + v[1] == v[1] + 2.0 - v[2]^2
		@test v[1]^2 + 8.81 == 9.81
	end
end

@testset "duals" begin
	@basis 4

	for dual in [flipdual, hodgedual, poincaredual]
		@test dual(v1 + 2v2) == dual(v1) + 2dual(v2)
	end

	for dim in 0:5, k in 0:dim
		a, b = KVector{dim,k}.(eachcol(rand(-5:5, ncomponents(dim, k), 2)))
		I = unit_pseudoscalar(dim)
		@test a ∧ hodgedual(b) == a ⊙ ~b * I
	end
end