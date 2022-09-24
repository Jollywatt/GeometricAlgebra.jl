using GeometricAlgebra:
	show_signature,
	show_multivector

@testset "show methods" begin
	io = IOBuffer()
	for a in [
		Blade{(1,1)}(0b00 => 42),
		Blade{(1,)}(0b1 => 1.1),
		Blade{(1,1,1)}(0b111 => true),
	], b in [a, Multivector(a), MixedMultivector(a)]
		@test isnothing(show(io, b))
		@test isnothing(show(io, MIME("text/plain"), b))
		@test repr(b) isa String
	end
end

@testset "pretty-printed metric signatures" begin
	sig = (+1,-1,-1,-1)
	prettysig = sprint(show_signature, sig)

	for T in [Blade, Multivector, MixedMultivector, HomogeneousMultivector]
		@test contains(sprint(show, T{sig}), prettysig)
		@test contains(sprint(show, Val{T{sig}}), prettysig)
		@test sprint(show, T) == string(nameof(T))
	end

	@test contains(sprint(show, Blade{sig,2}), prettysig)
	@test contains(sprint(show, Blade{sig,2,Int}), prettysig)
	@test contains(sprint(show, Blade{sig,K,Int} where {K}), prettysig)

	@test startswith(sprint(show, MIME("text/plain"), Blade{3}), "Blade{3")
	@test startswith(sprint(show, MIME("text/plain"), Blade{(1,1,1)}), "Blade{⟨+++⟩")
end

@testset "multivector printing" begin
	@basis 3

	@test sprint(show_multivector, 1000v1 + v2 + 0.001v3) == """
		1000.0   v1
		   1.0   v2
		   0.001 v3"""

	@test repr(zero(v1 + v2)) == "0"
	@test repr(zero(v1 + v12)) == "0"

	@test repr(1 + v1 + v2 + 3v12) == "1 + (1v1 + 1v2) + (3v12)"
end