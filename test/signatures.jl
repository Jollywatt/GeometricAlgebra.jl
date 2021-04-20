using GeometricAlgebra

@testset "indexing" begin
	sigs = [sig for s ∈ [(1,1,1), (x=1,y=1,z=1)] for sig ∈ [s, OffsetSignature(s,1:3)]]
	@testset "integer indices with $sig" for sig ∈ sigs
		x, y, z = basis(sig)
		@test x[1] == 1
		@test x[2] == 0
		@test (x + 4y)[2] == 4
		@test (1 + 2*x*y)[2,1] == -2
	end

	labelled3d = (x=1,y=1,z=1)
	@testset "symbol indices with $sig" for sig ∈ [labelled3d, OffsetSignature(labelled3d,11:13)]
		x, y, z = basis(sig)
		@test x[:x] == 1
		@test y[:x,:y] == 0
		@test (x + 4y)[:y] == 4
		@test (1 + 2*x*y)[:y,:x] == -2
	end

	@testset "offset indices" begin
		t, x, y, z = basis(OffsetSignature((t=-1, x=1, y=1, z=1),0:3))
		@test t[0] == x[1] == 1
		@test (x + 2y)[2] == 2
		@test (1 + 2*x*y)[2,1] == -2
		@test (t*x*y*z)[0,1,2,3] == 1
	end

end

@testset "@basis" begin
	@basis x y z=-1
	@test x isa Blade{signature(x),1}
	@test x^2 == 1
	@test z^2 == -1

	@basisall x y
	@test xy == -yx
end