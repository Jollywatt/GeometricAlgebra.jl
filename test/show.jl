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
@test isnothing(show(stdout, MIME("text/plain"), x))
@test isnothing(show(stdout, MIME("text/plain"), x + y))
@test isnothing(show(stdout, MIME("text/plain"), 1 + y))

@testset "pretty-printed types" begin
	for sig ∈ [(1, 1, 1), (x=1, y=1, z=1)]
		pretty_sig = GeometricAlgebra.show_signature(sig)
		@test sprint(show, Blade{sig,2,0b011,Float64}) == "Blade{$pretty_sig, 2, 0b011, Float64}"
		@test sprint(show, Blade{sig,2,:bits}) == "Blade{$pretty_sig, 2, :bits, T} where T"
		@test sprint(show, Blade{sig,2,0b101,Float64} where sig) == "Blade{sig, 2, 0b101, Float64} where sig"
		@test sprint(show, Multivector{sig,k,S} where {k,S}) == "Multivector{$pretty_sig, k, S} where {k, S}"
	end
	@test sprint(show, MixedMultivector) == "MixedMultivector"
end