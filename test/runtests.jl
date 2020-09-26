using Test
using GeometricAlgebra

ublade_types = [UInt8, UInt64, Vector{Int}]

@testset "lindex <-> ublade" begin
	@testset "$T" for T ∈ ublade_types
		@testset "$((; n, k))" for n ∈ 1:8, k ∈ 0:n
			for i ∈ 1:(n < 5 ? 1 : 11):binomial(n, k)
				@test GeometricAlgebra.ublade2lindex(GeometricAlgebra.lindex2ublade(T, k, i)) == i
			end
		end
	end
end

sigs = [
	(1,), (1,1), (1,1,1), (-1,-1,-1), (-1,1,1,1), (-1,-1,1,1), (0,1,1,1), (1,0,-1),
	(x=1,y=1,z=1), (t=1,x=-1,y=-1,z=-1)
]

@testset "ublade conversion" begin
	for a2b ∈ [
		0b1 => [1],
		0b10 => [2],
		0b11 => [1,2],
		UInt32(0b11001) => [1,4,5],
		UInt128(0b11101) => [1,3,4,5],
		# BitVector([1,1,0,1]) => [1,2,4], # havent implemented bv conversions yet
		# BitVector([1,1,0,1]) => 0b1011,
	], (a, b) ∈ [a2b, reverse(a2b)] # check both directions
		bfroma = GeometricAlgebra.convert_ublade(:sig, typeof(b), a)
		@test typeof(bfroma) == typeof(b)
		@test bfroma == b
	end
end

import GeometricAlgebra: ubladeprod
@testset "ubladeprod" begin
	@test ubladeprod((1,1,1), [1,2,1,2]) == (-1, [])
	@test ubladeprod((1,1,1), [1,2], [1,2]) == (-1, [])
	@test ubladeprod((1,1,1), [1,2], [3]) == (1, [1,2,3])
	@test ubladeprod((1,1,1), 0b11, 0b11) == (-1, 0b0)
	@test ubladeprod((-1,-1,-1), [1,3,1]) == (1, [3])
	@test ubladeprod((1,0,1), [2], [1,2]) == (0, [1])
end

# @testset "basis" begin
# 	sig = (1,1,1)
# 	x = Blade{sig}(1, 0b001)
#	# doesn't work yet because == requires implicit convert_ublade to compare
# 	@test basis(Blade{sig, Int, UInt8, 1}, 1) == x
# 	@test basis(Blade{sig, Float64, Vector{Int}, 2}, 1) == Blade{sig}(1, 0b011)
# end


@testset "Hodge star operator" begin
	for sig ∈ sigs
		n = dim(sig)
		s = prod(sig) # parity of signature, det(g)
		for η ∈ GeometricAlgebra.fullbasis(sig)
			k = grade(η)
			@test ★(★(η)) == (-1)^(k*(n - k))*s*η
			# identity: https://ncatlab.org/nlab/show/Hodge+star+operator#BasicProperties
		end
	end
end