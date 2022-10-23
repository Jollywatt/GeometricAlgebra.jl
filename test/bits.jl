using GeometricAlgebra:
	bits_to_indices,
	indices_to_bits,
	bits_of_grade,
	componentbits,
	componentindex,
	sign_from_swaps,
	factor_from_squares,
	geometric_prod_bits


@testset "bits <-> indices" begin
	for bits ∈ [0b000, 0b101, 0b110, 0b1101001, rand(UInt16, 3)..., rand(UInt32, 3)...]
		@test indices_to_bits(bits_to_indices(bits)) == bits
	end
end

@testset "bit permutations" begin
	@test first(GeometricAlgebra.BitPermutations{UInt8}(3)) === UInt8(0b111)
	@test first(GeometricAlgebra.BitPermutations{UInt32}(3)) === UInt32(0b111)
	@test collect(Iterators.take(bits_of_grade(2), 3)) == [0b011, 0b101, 0b110]
	@test collect(Iterators.take(bits_of_grade(5), 3)) == [0b011111, 0b101111, 0b110111]
	@test Iterators.take(bits_of_grade(10), 10) |> collect |> last == 0b11111111101
end

@testset "bits <-> kvector index" begin
	cases = [
		(3, 0, 1) => 0b000,
		(3, 1, 1) => 0b001,
		(3, 1, 2) => 0b010,
		(3, 2, 1) => 0b011,
		(3, 1, 3) => 0b100,
		(3, 2, 2) => 0b101,
		(3, 2, 3) => 0b110,
		(3, 3, 1) => 0b111,
		(8, 7, 1) => 0b01_111_111,
		(8, 7, 2) => 0b10_111_111,
	]
	@testset "no. $(i) of grade $k: $bits" for ((n, k, i), bits) ∈ cases
		@test componentbits(Val(n), Val(k))[i] == bits
		@test componentindex(KVector{n,k}, bits) == i
	end
end

@testset "bits <-> multivector index" begin
	for dim ∈ 1:2:8, i ∈ 1:2:2^dim
		@test componentindex(Multivector{dim}, componentbits(Val(dim))[i]) == i
	end
end

@testset "multivector slices" begin
	for n ∈ 0:8, k ∈ 0:n
		slice = componentindex(Multivector{n}, KVector{n,k})
		@test componentbits(Val(n))[slice] == componentbits(Val(n), Val(k))
	end
end

@testset "sign_from_swaps()" begin
	# note that typographically bits are written 'backwards'; v124 == 0b1011
	@test sign_from_swaps(0b01, 0b10) == +1
	@test sign_from_swaps(0b10, 0b01) == -1
	@test sign_from_swaps(0b11, 0b01) == -1
	@test sign_from_swaps(0b111, 0b01) == +1
end

@testset "factor_from_squares()" begin
	@test factor_from_squares((+1, +1, +1), 0b000) == +1
	@test factor_from_squares((+1, +1, +1), 0b111) == +1
	@test factor_from_squares((-1, +1, +1), 0b001) == -1
	@test factor_from_squares((-1, +1, +1), 0b010) == +1
	@test factor_from_squares((-1, +1, +1), 0b111) == -1
	@test factor_from_squares((-1, -1, -1), 0b101) == +1
	@test factor_from_squares((-1, -1, -1), 0b111) == -1
end

@testset "geometric_prod_bits()" begin
	@test geometric_prod_bits((1, 1, 1), 0b001, 0b010) == (+1, 0b011)
	@test geometric_prod_bits((1, 1, 1), 0b010, 0b001) == (-1, 0b011)
	@test geometric_prod_bits((1, 1, 1), 0b011, 0b001) == (-1, 0b010)
	@test geometric_prod_bits((0, 1, 1), 0b001, 0b101) == (0, 0b100)
end