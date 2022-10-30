using GeometricAlgebra:
	bits_of_grade

@testset "scalar *" begin
	a = BasisBlade{(1,1)}(0b01 => 10)

	for b in [a, Multivector(a)]
		@test 3a == a*3.0
		@test a/10 == 10\a
		@test -a == a/(-1.0)
	end

	@test a//5 === BasisBlade{(1,1)}(0b01 => 2//1)
end

@testset "+" begin
	v = basis(3)
	bi = basis(3, grade=2)

	@test v[1] + v[2] isa Multivector{3,1}
	@test bi[1] + 2.5bi[2] isa Multivector{3,2}
	@test v[1] + bi[2] isa Multivector{3,0:3}
	@test sum(Multivector.(bi)) isa Multivector{3,2}
	@test v[1] - v[2] == -v[2] + v[1]

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

	@test geometric_prod(v[1], 2) == 2v[1]

	Ts = [BasisBlade, Multivector]
	for T1 in Ts, T2 in Ts
		@test v[1]v[2] == T1(v[1])T2(v[2]) == T2(v[1])T1(v[2])
	end
end

@testset "scalar product" begin
	v = basis((-1,+1,+1,+1))

	@test v[1] ⊙ v[1] == -1
	@test v[1] ⊙ v[2] == 0
	@test v[1] ⊙ (1 + v[1]v[2]) isa Real

	@test scalar_prod(v[2] + v[3], v[2] + v[3]) == 2
	@test scalar_prod(v[1] + v[2], v[3]) == 0
	@test scalar_prod(10 + 2v[3], 20 + v[3]) == 202
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
	@basisall 3

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

	for dual in [flipdual, hodgedual, poincaredual]
		@test dual(5v[1]) == 5dual(v[1])
		@test dual(v[1] + 2v[2]) == dual(v[1]) + 2dual(v[2])
		@test dual(v[1] + 10v[1]v[2]v[3]) == dual(v[1]) + 10dual(v[1]v[2]v[3])
	end

	for dim in 0:5, k in 0:dim
		a, b = Multivector{dim,k}.(eachcol(rand(-5:5, ncomponents(dim, k), 2)))
		I = flipdual(one(a)) # unit pseudoscalar
		@test a ∧ hodgedual(b) == a ⊙ ~b * I

		m = Multivector{dim,0:dim}(rand(-5:5, ncomponents(dim)))
		@test hodgedual(m) == ~m*I
		@test flipdual(flipdual(m)) == m
	end
end