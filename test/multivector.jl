@testset "promotion" begin
	sig = (x=1,y=1,z=1)
	xs = Any[
		Complex(42)
		Blade{sig}(42, 0b101)
		Blade{sig}(42.0, [1, 3])
		Blade{sig}(42//1, [:x, :z])
		Multivector{sig,0,Vector{Float32}}([42])
		Multivector{sig,0,Vector{Int16}}([42])
		Multivector{sig,2,Dict{UInt,Int8}}(Dict(0b101 => 42))
		Multivector{sig,2,Dict{Vector{Symbol},Int}}(Dict([:x, :z] => 42))
		convert(MixedMultivector{sig,Dict{Vector{Symbol},Int128}}, 42)
		convert(MixedMultivector{sig,Vector{Int32}}, 42)
	]
	for a ∈ xs, b ∈ xs
		T = promote_type(typeof.((a, b))...)
		@test promote(a, b) isa Tuple{T,T}
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


@testset "component access" begin
	t, x, y, z = basis((t=-1, x=1, y=1, z=1))
	a = 1 + 2t + 3*x*z - x*y*z

	@testset "integer slots" begin
		@test a[] == 1
		@test a[1] == 2
		@test a[2,4] == 3
		@test a[2] == 0
		@test a[4,2] == -3
		@test a[2,3,4] == a[4,2,3] == -1
	end

	@testset "symbol slots" begin
	    @test a[:t] == 2
	    @test a[:x] == 0
	    @test a[:z,:x] == -3
	end


end