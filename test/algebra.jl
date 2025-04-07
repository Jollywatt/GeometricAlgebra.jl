using GeometricAlgebra:
	bits_of_grade,
	add!

using StaticArrays

@testset "==" begin
	v = basis(2)
	@test v[1] + 0v[2] == v[1]
	@test 0v[1] == 0
	@test 0v[1] + 0 == 0
	@test iszero(1 + v[2] - v[2] - 1)
	@test isone(v[1] + -v[1] + 1)
	@test !isone(v[1])
	@test !isone(Multivector(v[1]))
end

@testset "≈" begin
	v = basis(3)
	@test v[1] + √eps(1.0) ≈ v[1] rtol=1e-6
	@test 1e-10 + v[1] ≈ v[1] rtol=1e-10
	@test sum(v)/1e16 ≈ 0 atol=1e-10 
	@test sin(π*v[1]) ≈ 0 atol=1e-10
end

@testset "scalar *" begin
	a = BasisBlade{(1,1)}(10, 0b01)

	for b in [a, Multivector(a)]
		@test 3a == a*3.0
		@test a/10 == 10\a
		@test -a == a/(-1.0)
	end

	@test a//5 === BasisBlade{(1,1)}(2//1, 0b01)
end

@testset "add!" begin
	u = Multivector{3,1}([1,2,3])
	add!(u, 100, 0b001)
	@test u == Multivector{3,1}([101,2,3])

	v = copy(u)
	add!(u, 100, 0b111)
	@test u == v

	add!(u, u)
	@test u == 2v

	# non-mutable component type:
	# add! should return a new instance of identical type
	u = Multivector{3,1}(SA[1,2,3])
	v = add!(u, 100, 0b001)
	@test u !== v
	@test typeof(u) === typeof(v)

	m = Multivector{3,0:3}(zeros(Int, 8))
	add!(m, u)
	@test m == u
end

@testset "+" begin
	v = basis(3)
	bi = basis(3, 2)

	@test v[1] + v[2] isa Multivector{3,1}
	@test bi[1] + 2.5bi[2] isa Multivector{3,2}
	@test v[1] + bi[2] isa Multivector{3,0:3}
	@test sum(Multivector.(bi)) isa Multivector{3,2}
	@test v[1] - v[2] == -v[2] + v[1]
	@test (v[1] + v[2]) - (v[2] - v[3]) == v[1] + v[3]

	@testset "with scalars" begin
		@test v[1] + 7.5 == 7.5 + v[1]
		@test (v[1] + v[2]) + 7 == v[1] + (v[2] + 7)
		@test (bi[1] + 8) - 8 == bi[1]
		@test 1 - v[1] == -(v[1] - 1)
		@test Multivector(v[1]) + 0.5 == v[1] + 1//2
	end

	@test (1:3)'v == 1v[1] + 2v[2] + 3v[3]
end

@testset "*" begin
	v = basis((-1,+1,+1,+1))

	@test v[1]v[2] == -v[2]v[1]
	@test v[1]v[1] == -1
	@test (v[1] + 7v[3])v[2] == v[1]v[2] + 7v[3]v[2]
	@test (v[1]v[1] + v[2])v[2] == v[2]v[2] - v[2]
	@test v[1]v[2]v[3] == v[2]v[3]v[1] == v[3]v[1]v[2]
	@test 4 == 2*(v[2] + v[3])*(v[2] + v[3])

	@test geometric_prod(v[1], 2) === 2v[1] === geometric_prod(2, v[1])
	@test geometric_prod(6, 7) === 42

	Ts = [BasisBlade, Multivector]
	for T1 in Ts, T2 in Ts
		@test v[1]v[2] == T1(v[1])T2(v[2]) == T2(v[1])T1(v[2])
	end
end

@testset "scalar product" begin
	v = basis((-1,+1,+1,+1))

	@test v[1] ⊙ v[1] === -1
	@test v[1] ⊙ v[2] === 0
	@test v[1] ⊙ (v[1]v[2]) === 0
	@test v[1] ⊙ (1 + v[1]v[2]) isa Real

	@test v[1] ⊙ 7 === 0
	@test 6⊙7 === 42

	@test scalar_prod(v[2] + v[3], v[2] + v[3]) === 2
	@test scalar_prod(v[1] + v[2], v[3]) === 0
	@test scalar_prod(10 + 2v[3], 20 + v[3]) === 202
end

@testset "∧" begin
	v = basis(3)

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
	@basis 3 allperms=true

	@test v123⨽v32 == v1 == ~(v23⨼v321)
	@test iszero(7⨽(v1 + v2))
	@test v2⨽(10 + v123) == 10v2
	@test 2⨼2 == 4
end

@testset "^" begin
	v = basis((-1,+1,+1,+1))
	
	@test v[1]^2 + 1 == 0
	@test v[1]^4003 == -v[1]

	@test (v[2] + v[3])^2 == 2
	@test (1 + v[1])^2 == 2v[1]

	a = 1 + 2v[2] + 3v[1]v[2]
	@test a^30 == ((a^2)^3)^5

	@test (v[1]v[2] + v[2]v[3])^10 == ((v[1]v[2])^2 + (v[2]v[3])^2)^5
end

@testset "reversion" begin
	v = basis((+1,-1,-1,0))
	
	@test reversion(v[1]) == v[1]
	@test reversion(v[1]v[2]) == v[2]v[1]
	@test reversion(3v[1]v[2] + v[3]v[4]) == 3v[2]v[1] + v[4]v[3]
	@test reversion(1 + 2v[1] + 3v[2]v[3] + 4v[1]v[2]v[3]) == 1 + 2v[1] + 3v[3]v[2] + 4v[3]v[2]v[1]

end

@testset "duals" begin
	v = basis(4)

	for dual in [flipdual, rdual, ldual, hodgedual, invhodgedual]
		@test dual(5v[1]) == 5dual(v[1])
		@test dual(v[1] + 2v[2]) == dual(v[1]) + 2dual(v[2])
		@test dual(v[1] + 10v[1]v[2]v[3]) == dual(v[1]) + 10dual(v[1]v[2]v[3])
	end

	for dim in 0:5, k in 1:dim
		a, b = Multivector{dim,k}.(eachcol(rand(-5:5, ncomponents(dim, k), 2)))
		I = flipdual(one(a)) # unit pseudoscalar
		@test a ∧ hodgedual(b) == a ⊙ ~b * I

		m = Multivector{dim,0:dim}(rand(-5:5, ncomponents(dim)))
		@test hodgedual(m) == ~m*I
		@test flipdual(flipdual(m)) == m
	end

end

@testset "hodgedual, invhodgedual" begin
	for sig in [2, 3, "++++", "-+++", "--++", 5]
		V = basis(sig, :all)
		@test hodgedual.(invhodgedual.(V)) == V
		@test invhodgedual.(hodgedual.(V)) == V
	end

	for sig in ["+++0", "++0-"]
		V = basis(sig, :all)
		@test all(@. ifelse(iszero(hodgedual(V)), hodgedual(invhodgedual(V)) == V, invhodgedual(hodgedual(V)) == V))
	end
end

@testset "∨" begin
	for sig in ["+++", "++-", "+00", "---", "-+++", "+--0"]
		V = basis(sig, :all)
		for dual ∈ [ldual, rdual]
			@test all(dual(a∨b) == dual(a)∧dual(b) for a in V, b in V)
			@test all(dual(a)∨dual(b) == dual(a∧b) for a in V, b in V)
		end
	end
end

@testset "sandwich_prod" begin
	@basis 3

	R = exp(π/4*v12)
	@test sandwich_prod(R, v2) ≈ v1
	@test grade(sandwich_prod(R, v2)) == 1
	@test sandwich_prod(R, 7) ≈ 7
	@test sandwich_prod(2, R) == 4R
	@test sandwich_prod(4, 2) == 4*2*4
end

@testset "outermorphism" begin
	@basis 3 scalar=true

	𝟙 = [1 0 0; 0 1 0; 0 0 1]
	@test outermorphism(𝟙, 7) == 7
	@test outermorphism(𝟙, 7v) == 7v
	@test outermorphism(𝟙, 7v1) == 7v1
	@test outermorphism(𝟙, 7v1 + v23) == 7v1 + v23

	A = [
	1 0 0
	0 2 1
	0 0 1
	]
	@test outermorphism(A, v1) == v1
	@test outermorphism(A, v2) == 2v2
	@test outermorphism(A, v3) == v2 + v3
	@test outermorphism(A, v23) == 2v2∧(v2 + v3)
	@test outermorphism(A, v123) == 2v123
end

@static if VERSION >= v"1.8"
	@testset "inferability" begin
		@basis 3 scalar=true

		a, b = v1 + v2, v12 + v23

		@inferred v1 + v12
		@inferred v + v123
		@inferred a + b

		@inferred a*a
		@inferred b*b
		@inferred a*b
		@inferred Multivector(v)*Multivector(v123)
		@inferred Multivector(v1)*Multivector(v23)
		@inferred Multivector(v12)*Multivector(v23)

		@inferred (v + v12) + v23
		@inferred (v + v12) + v3
		@inferred (v + v12)*v3

		@inferred ~(rdual(a)*b + v13) + b

		@inferred v12∧v3
		@inferred a∧b
		@inferred (a⨼b)∧(v12⋅b)*I
	end
end

@testset "embed" begin
	a = Multivector{3,1}(1, 2, 3)
	@test embed(4, a) == Multivector{4,1}(1, 2, 3, 0)
	@test embed(2, a) == Multivector{2,1}(1, 2)

	a = randn(Multivector{3,0:3})
	@test embed(3, embed(Cl(4,1), a)) == a
end
