using GeometricAlgebra:
	via_matrix_repr,
	inv_matrix_method

@testset "inverses" begin
	v = basis((-1,+1,+1,+1))
	@test inv(v[1])v[1] == 1
	@test inv(5v[1]v[2]) == inv(v[2])inv(v[1])inv(5)

	u = sum(v)
	@test 1/u == u\1 == u/scalar(u^2) == scalar(u^2)\u

	@testset "dimension $dim" for dim in 1:5
		mixed_sig = (-1,+1,+1,+1,0)
		for sig in [
			dim,
			mixed_sig[1:dim],
		], trials in 1:5
			a = Multivector{sig,0:dim}(rand(2^dim))
			@test a/a ≈ 1 rtol=1e-5
			@test 1 ≈ a\a rtol=1e-5
			@test inv_matrix_method(inv(a)) ≈ a rtol=1e-6
		end
	end
end

@testset "sqrt" begin
	@basis 3 scalar=true

	# sqrt of scalar
	@test sqrt(4v) === 2.0v
	@test sqrt(4Multivector(v)) == 2
	@test sqrt(-v) == I

	@test sqrt(I)^2 ≈ I

	u = Multivector{3,1}([1,2,3]/4)
	@test sqrt(u)^2 ≈ u
	@test sqrt(u*I)^2 ≈ u*I
	@test sqrt(u + I)^2 ≈ u + I rtol=sqrt(eps())
end

@testset "exp" begin
	v = basis((1,1))
	@test exp(10000*2pi*v[1]v[2]) ≈ 1

	@test exp(0v[1]) == 1
	@test exp(3v[1]) ≈ cosh(3) + sinh(3)v[1]
	@test exp(3v[1]v[2]) ≈ cos(3) + sin(3)v[1]v[2]

	for dim in 1:5, trials in 1:5
		a = Multivector{dim,0:dim}(big.(randn(2^dim)))
		@test exp(a)exp(-a) ≈ 1
		@test inv(exp(a)) ≈ exp(-a)
	end

	@test exp((3 + 0im)v[1]) ≈ exp(3v[1])
	@test exp((3 + 2im)v[1]) ≈ GeometricAlgebra.exp_series(Multivector((3 + 2im)v[1]))
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
		for a in [5v[1]v[2], v[1]v[2] + v[2]v[3], rand(length(v))'v, im*v[1]]
			@test exp(a) ≈ cosh(a) + sinh(a)
			@test tan(a) ≈ sin(a)/cos(a)
			@test sin(a) ≈ sin(asin(sin(a)))
		end
	end
end

@testset "compact representations" begin
	# if grades fit in a smaller subalgebra,
	# restrict result to that subalgebra

	for a in [
		Multivector{4,0}([7])
		Multivector{Cl(1,3),0:4:4}([4, 5])
		Multivector{4,0:2:4}(fill(1, 8))
		Multivector{4,0:4}(fill(1, 16))
	]

		@test grade(sqrt(a)) === grade(a)
		@test grade(exp(a)) === grade(a)
		@test grade(log(a)) === grade(a)
		@test grade(sin(a)) === grade(a)
	end

	m02 = Multivector{4,(0,2)}(1:7)

	@test sqrt(m02) |> grade === 0:2:4
	@test exp(m02) |> grade === 0:2:4

	@test via_matrix_repr(sin, Multivector{3,0:2:2}(1:4)) |> grade === 0:2:2
end


@testset "blade factorisation" begin
	@testset "blades" begin
		@testset for T in (Float64, ComplexF64)
			@testset "sig: $sig" for sig in 1:5
				@testset "grade: $k" for k in 1:dimension(sig)
					@testset for b in basis(sig, k)
						b *= randn(T)
						factors = factorblade(b)
						@test wedge(factors...) ≈ b
					end
				end
			end
		end
	end

	@testset "multivectors" begin
		@testset for T in (Float64, ComplexF64)
			@testset "sig: $sig" for sig in [2, 3, 4, 5, Cl(1,1), Cl(1,3)]
				@testset "grade: $k" for k in 1:dimension(sig)
					for _ in 1:100
						blade = wedge(randn(Multivector{sig,1,Vector{T}}, k)...)
						factors = factorblade(blade)
						@test wedge(factors...) ≈ blade
					end
				end
			end
		end
	end
end
