using GeometricAlgebra:
	show_signature,
	show_multivector
using StaticArrays

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

	for sig in [3, (1,1,1), Cl(2), Cl(3,0,1)]
		@test startswith(sprint(show, MIME("text/plain"), BasisBlade{sig}), "BasisBlade{"*sprint(show_signature, sig))
	end
end

@testset "multivector printing" begin
	v = basis(3)

	@test sprint(1000v[1] + v[2] + 0.001v[3]) do io, a
		show_multivector(io, a, inline=false)
	end == """
		1000.0   v1
		   1.0   v2
		   0.001 v3"""

	@test sprint(1000v[1] + v[2] + 0.001v[3]) do io, a
		show_multivector(io, a, inline=true, compact=false)
	end == "1000.0 v1 + 1.0 v2 + 0.001 v3"

	@test sprint(1000v[1] + v[2] + 0.001v[3]) do io, a
		show_multivector(io, a, inline=true, compact=true)
	end == "1000.0v1 + v2 + 0.001v3"

	@test sprint(1 + 1v[1] + 2v[2] + v[1]v[2]v[3]) do io, a
		show_multivector(io, a, inline=false, groupgrades=true)
	end == """
		1
		1 v1 + 2 v2
		1 v123"""

	@test sprint(1 + 1v[1] + 2v[2] + v[1]v[2]v[3]) do io, a
		show_multivector(io, a, inline=false, groupgrades=false)
	end == """
		1 v
		1 v1
		2 v2
		1 v123"""

	@test sprint(1 + 1v[1] + 2v[2] + v[1]v[2]v[3]) do io, a
		show_multivector(io, a, inline=true, groupgrades=false)
	end == "1 + 1 v1 + 2 v2 + 1 v123"

end

@testset "repr parseability" begin
	v = basis(3)
	reflect(eq, x) = eq(eval(Meta.parse(repr(x))), x)

	@test reflect(===, v[1])
	@test reflect(===, 5.5v[1])
	@test reflect(===, 4//v[1])

	@test reflect(==, Multivector{3,1}([1,2,3]))
	@test reflect(==, Multivector{3,2}(SVector(1,2,3)))
	@test reflect(==, Multivector{3,3}(MVector(42)))
end