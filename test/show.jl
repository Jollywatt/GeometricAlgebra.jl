using GeometricAlgebra:
	show_signature

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
end