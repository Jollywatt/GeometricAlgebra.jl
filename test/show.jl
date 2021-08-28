using GeometricAlgebra

v1, v2 = basis((1, 1))
v12 = v1*v2

x, y = basis((x=1, y=1))
xy = x*y

@testset "repr parseability" begin
	isparseable(a) = eval(Meta.parse(repr(a))) == a
	
	@test isparseable(4v1)
	@test isparseable(0v1)
	@test isparseable((1 + im)*v1)
	@test isparseable(v1 - 2v2)
	@test isparseable(4 + v12)

	@test isparseable(10x)
	@test isparseable(0*x)
	@test isparseable(0.0x)
	@test isparseable(1 + 2x + 3xy)
end

# just test against errors
for to_show ∈ (x, x + y, 1 + y, Any[100x, 2x*y, 1 + x, 0 + 0y], typeof(x), typeof(x + y))
	@test isnothing(show(stdout, MIME("text/plain"), to_show))
end

@testset "pretty-printed types" begin
	for sig ∈ [(1, 1, 1), (x=1, y=1, z=1)]
		pretty_sig = GeometricAlgebra.show_signature(sig)
		compare(a, b) = startswith(replace(a, r"\s+"=>""), replace(b, r"\s+"=>""))
		@test compare(sprint(show, Blade{sig,2,Float64}), "Blade{$pretty_sig, 2, Float64}")
		@test compare(sprint(show, Blade{sig,2,Float64} where sig), "Blade{sig, 2, Float64} where sig")
		@test compare(sprint(show, Multivector{sig,k,S} where {k,S}), "Multivector{$pretty_sig")
	end
	@test sprint(show, MixedMultivector) == "MixedMultivector"
end

@testset "blade printing" begin
	@test repr(Blade{(1,1,1)}(42, 0b11)) == "42v12"
	@test repr(Blade{(x=1,y=1,z=1)}(42, 0b11)) == "42xy"
	@test repr(Blade{OffsetSignature((-1,1,1,1),0:3)}(42, 0b11)) == "42v01"
	@test repr(Blade{OffsetSignature((t=-1,x=1,y=1,z=1),0:3)}(42, 0b11)) == "42tx"
end