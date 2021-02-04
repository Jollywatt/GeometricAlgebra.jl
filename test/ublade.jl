using GeometricAlgebra:
	ublade2lindex, lindex2ublade,
	convert_ublade, ubladeprod

@testset "lindex ↔︎ ublade" begin
	@testset "ublade type $T" for T ∈ [UInt8, UInt, Vector{Int}]
		for n ∈ 1:8, k ∈ 0:n, i ∈ 1:(n < 5 ? 1 : 13):binomial(n, k)
			@test ublade2lindex(lindex2ublade(T, k, i)) == i
		end
	end
end

@testset "ublade conversion" begin
	sig = (t=-1, x=1, y=1, z=1)
	for a2b ∈ [
		0b0 => Int[],
		0b1 => [1],
		0b10 => [2],
		0b11 => [1,2],
		UInt32(0b1001) => [1,4],
		UInt128(0b1101) => [1,3,4],
		Int[] => Symbol[],
		[1,2,4] => [:t,:x,:z],
		0b0 => Symbol[],
		0b1110 => [:x,:y,:z],
	], (a, b) ∈ [a2b, reverse(a2b)] # check both directions
		bfroma = convert_ublade(sig, typeof(b), a)
		@test typeof(bfroma) == typeof(b)
		@test bfroma == b
	end
end


@testset "ubladeprod" begin
	@test ubladeprod((1,1,1), [1,2,1,2]) == (-1, [])
	@test ubladeprod((1,1,1), [1,2], [1,2]) == (-1, [])
	@test ubladeprod((1,1,1), [1,2], [3]) == (1, [1,2,3])
	@test ubladeprod((1,1,1), 0b11, 0b11) == (-1, 0b0)
	@test ubladeprod((-1,-1,-1), [1,3,1]) == (1, [3])
	@test ubladeprod((1,0,1), [2], [1,2]) == (0, [1])
	@test ubladeprod((x=1,y=1,z=1), [:x,:y], [:x,:y]) == (-1, [])
end