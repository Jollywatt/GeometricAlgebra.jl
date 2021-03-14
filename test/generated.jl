@testset "generated products" begin
	using GeometricAlgebra: geometric_prod, geometric_prod_gen, homogeneous_prod, homogeneous_prod_gen
	for sig ∈ [(1,1), (t=-1,x=+1,y=+1,z=+1), (+1,0,-1)], sample ∈ [-100:100, Float64]
		A, B = [MixedMultivector{sig}(rand(MersenneTwister(i), sample, 2^dimension(sig))) for i ∈ 1:2]
		U, V = [Multivector{sig,1}(rand(MersenneTwister(i), sample, dimension(sig))) for i ∈ 1:2]

		for a ∈ (A, U), b ∈ (B, V)
			@test geometric_prod_gen(a, b) == geometric_prod(a, b)
			for k ∈ 1:dimension(sig)
				@test homogeneous_prod_gen(a, b, Val(k)) == homogeneous_prod(a, b, k)
			end
		end
	end
end