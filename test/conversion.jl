let
	sig = (x=1,y=1,z=1)
	x, y, z = basis(Blade{sig,1,Int,UInt})

	@testset "complex, real, imag" begin
		for u ∈ (x, x*y + y*z, 1 + x*z)
			@test complex(typeof(u)) == typeof(complex(u))
			@test real(complex(typeof(u))) == typeof(u)
			@test real(complex(u)) == imag(u*im)
			z = u + 2im*u
			@test z == real(z) + imag(z)im
		end
	end

	@testset "float, big" begin
		for u ∈ (x, x*y + y*z, 1 + x*z)
			@test typeof(float(u)) == float(typeof(u))
			@test eltype(float(u)) == float(eltype(u))
			@test typeof(big(u)) == big(typeof(u))
			@test eltype(big(u)) == big(eltype(u))
		end
	end
end