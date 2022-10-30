using GeometricAlgebra:
	show_signature,
	show_multivector


@testset "pretty-printed metric signatures" begin
	sig = (+1,-1,-1,-1)
	prettysig = sprint(show_signature, sig)

	for T in [BasisBlade, Multivector]
		@test contains(sprint(show, T{sig}), prettysig)
		@test contains(sprint(show, Val{T{sig}}), prettysig)
		@test sprint(show, T) == string(nameof(T))
	end

	@test contains(sprint(show, BasisBlade{sig,2}), prettysig)
	@test contains(sprint(show, BasisBlade{sig,2,Int}), prettysig)
	@test contains(sprint(show, BasisBlade{sig,K,Int} where {K}), prettysig)

	@test startswith(sprint(show, MIME("text/plain"), BasisBlade{3}), "BasisBlade{3")
	@test startswith(sprint(show, MIME("text/plain"), BasisBlade{(1,1,1)}), "BasisBlade{⟨+++⟩")
end

@testset "multivector printing" begin
	@basis 3

	@test sprint(1000v1 + v2 + 0.001v3) do io, a
		show_multivector(io, a, inline=false)
	end == """
		1000.0   v1
		   1.0   v2
		   0.001 v3"""

	@test sprint(1000v1 + v2 + 0.001v3) do io, a
		show_multivector(io, a, inline=true, compact=false)
	end == "1000.0 v1 + 1.0 v2 + 0.001 v3"

	@test sprint(1000v1 + v2 + 0.001v3) do io, a
		show_multivector(io, a, inline=true, compact=true)
	end == "1000.0v1 + v2 + 0.001v3"

	@test sprint(1 + 1v1 + 2v2 + v123) do io, a
		show_multivector(io, a, inline=false, groupgrades=true)
	end == """
		1
		1 v1 + 2 v2
		1 v123"""

	@test sprint(1 + 1v1 + 2v2 + v123) do io, a
		show_multivector(io, a, inline=false, groupgrades=false)
	end == """
		1 v
		1 v1
		2 v2
		1 v123"""

	@test sprint(1 + 1v1 + 2v2 + v123) do io, a
		show_multivector(io, a, inline=true, groupgrades=false)
	end == "1 + 1 v1 + 2 v2 + 1 v123"

	@test repr(zero(v1 + v2)) == "0"
	@test repr(zero(v1 + v12)) == "0"

	@test repr(1 + v1 + v2 + 3v12) == "(1) + (v1 + v2) + (3v12)"
end