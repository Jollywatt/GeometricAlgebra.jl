using Test, Random
using GeometricAlgebra

#=
TODO: test indexing
=#

ublade_types = [UInt8, UInt64, Vector{Int}]

@testset "lindex/ublade conversion" begin
	@testset "$T" for T ∈ ublade_types
		@testset "$((; n, k))" for n ∈ 1:8, k ∈ 0:n
			for i ∈ 1:(n < 5 ? 1 : 11):binomial(n, k)
				@test GeometricAlgebra.ublade2lindex(GeometricAlgebra.lindex2ublade(T, k, i)) == i
			end
		end
	end
end

sigs_tuple = [(1,), (1,1), (1,1,1), (-1,-1,-1), (-1,1,1,1), (1,0,-1), (0,)]
sigs_namedtuple = [(x=1,y=1,z=1), (t=-1,x=1,y=1,z=1), (ε=0,)]

sigs = [
	(1,), (1,1,1), (-1,-1,-1), (-1,1,1,1), (-1,-1,1,1), (1,0,-1),
	(x=1,y=1,z=1), (t=1,x=-1,y=-1,z=-1)
]

@testset "convert_ublade" begin
	for a2b ∈ [
		0b1 => [1],
		0b10 => [2],
		0b11 => [1,2],
		UInt32(0b11001) => [1,4,5],
		UInt128(0b11101) => [1,3,4,5],
	], (a, b) ∈ [a2b, reverse(a2b)] # check both directions
		bfroma = GeometricAlgebra.convert_ublade(:sig, typeof(b), a)
		@test typeof(bfroma) == typeof(b)
		@test bfroma == b
	end

	sig = (a=1, b=1, c=0)
	for a2b ∈ Pair{Vector{Int},Vector{Symbol}}[
		[] => [], [1] => [:a], [2, 3] => [:b, :c]
	], (a, b) ∈ [a2b, reverse(a2b)]
		bfroma = GeometricAlgebra.convert_ublade(sig, typeof(b), a)
		@test typeof(bfroma) == typeof(b)
		@test bfroma == b
	end
end

import GeometricAlgebra: ubladeprod
@testset "ubladeprod" begin
	@test ubladeprod((1,1,1), [1,2,1,2]) == (-1, [])
	@test ubladeprod((1,1,1), [1,2], [1,2]) == (-1, [])
	@test ubladeprod((1,1,1), [1,2], [3]) == (1, [1,2,3])
	@test ubladeprod((1,1,1), 0b11, 0b11) == (-1, 0b0)
	@test ubladeprod((-1,-1,-1), [1,3,1]) == (1, [3])
	@test ubladeprod((1,0,1), [2], [1,2]) == (0, [1])
	@test ubladeprod((x=1,y=1,z=1), [:x,:y], [:x,:y]) == (-1, [])
end

@testset "basis() and equality" begin
	sig = (1,1,1)
	x = Blade{sig}(1, 0b001)
	@test basis(Blade{sig, Int, UInt8, 1}, 1) == x
	@test basis(Blade{sig, Float64, Vector{Int}, 2}, 1) == Blade{sig}(1, 0b011)
	@test basis(sig, 1) == x
	@test basis(sig)[1] == x
end

@testset "equality" begin
	x, y = basis((1, 1))
	x_alts = [x, x + y - y, x + 0]
	for x1 ∈ x_alts, x2 ∈ x_alts
		@test x1 == x2
	end
	@test one(x) == 1
end

@testset "addn. & scalar mult." begin
	@testset "$sig" for sig ∈ sigs
		x = basis(sig, 1)
		@test x + 0 == x == 0 + x == 1 + x - 1
		@test x + x == 2x
		@test -x == x - 2x
		@test x/2 == 0.5x == 1//2*x == 2\x
	end
end


@testset "geometric prod" begin
	t, x, y, z = basis((-1,1,1,1))
	@testset "blades" begin
		@test x*y == -y*x
		@test x*x == 1
		@test t*t == -1
		@test (x*y)*(x*y) == x*y*x*y == -x*x*y*y == -1
		@test x != y
		@test x != x*y
		@test x != 0
	end
	@testset "(mixed) multivectors" begin
		@test (x + y)*z == x*z + y*z
		@test (x + 1)*x == x + 1
		@test (x - 2.)*y == x*y - y*2.
		@test (x + 1)*(x + 1) == 2x + 2
		@test x != x + 1
	end
end

@testset "other products" begin
	t, x, y, z = basis((t=-1,x=1,y=1,z=1))
	@test (x + y)∧z == x*z + y*z
	@test (x + y)∧x == y*x
	@test x∧(y + x) isa Multivector
	@test iszero(x∧x)

	sample = -100:100
	for sig ∈ sigs
		A, B, C = [MixedMultivector{sig}(rand(MersenneTwister(i), sample, 2^dim(sig))) for i ∈ 1:3]
		V = Multivector{sig}(rand(MersenneTwister(0), sample, dim(sig)))

		I = vol(A)

		@test V*A == V⨼A + V∧A
		@test A*V == A⨽V + A∧V

		iszero(I^2) && continue

		@test A⨼B == (A∧(B/I))*I
		@test A⨽B == I*((I\A)∧B)
		@test grade((A∧B)*C, 0) == grade(A*(B⨼C), 0)
		@test grade(C*(A∧B), 0) == grade((C⨽A)*B, 0)
		@test A⨼(B⨼C) == (A∧B)⨼C
		@test A⨼(B⨽C) == (A⨼B)⨽C
		@test reversion(A⨼B) == reversion(B) ⨽ reversion(A)
		@test reversion((reversion(A)/I)∧(B⨼C)) == (reversion(C)⨽reversion(B))∧(reversion(I)\A)

	end
end

@testset "Hodge star operator" begin
	for sig ∈ sigs
		n = dim(sig)
		s = prod(sig) # parity of signature, det(g)
		for η ∈ GeometricAlgebra.fullbasis(sig)
			k = grade(η)
			@test ★(★(η)) == (-1)^(k*(n - k))*s*η
			# identity: https://ncatlab.org/nlab/show/Hodge+star+operator#BasicProperties
		end
	end
end

@testset "inversion" begin
	# non-zero Euclidean vectors are always invertible
	for n ∈ 1:5
		v = Multivector{EuclideanSignature(n)}(rand(n) .+ 1e-3)
		@test v*inv(v) ≈ 1 ≈ inv(v)*v

		v = Multivector{EuclideanSignature(n)}(rand(1:10, n).//rand(1:5, n))
		@test isone(v*inv(v))
	end

	x, y, z = basis((1, 1, 1))

	# inverses in ℂ
	for (Re, Im) ∈ eachcol(rand(2, 10))
		c = complex(Re, Im)
		u = Re + Im*x*y
		ci = inv(c)
		@test ci.re + ci.im*x*y ≈ inv(u)
		@test inv(u)*u ≈ 1
	end

	i = basis(Blade{(-1,), Rational{Int}, UInt, 1}, 1)
	@test inv(2 + i)*(2 + i) == 1
end



@testset "trig funcs" begin
	t, x, y = basis((-1, 1, 1))

	for λ ∈ [0, 1, -0.5, π, -13, 11//5]
		@test exp(λ*x*y) == cos(λ) + x*y*sin(λ)
		@test exp(λ*x) == cosh(λ) + x*sinh(λ)
		for u ∈ [λ*x, λ*t, λ*x*y, λ*x*t, λ*x*y*t]
			@test tan(u) ≈ sin(u)/cos(u)
			@test tanh(u) ≈ sinh(u)/cosh(u)

			if abs(λ) <= 1
				@test tanh(atanh(u)) ≈ u
			end
		end
	end
end