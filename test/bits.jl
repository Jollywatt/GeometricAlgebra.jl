using GeometricAlgebra:
	bits_to_indices,
	indices_to_bits,
	BitPermutations,
	mv_index_to_bits,
	bits_to_mv_index,
	bits_to_mmv_index,
	mmv_index_to_bits,
	sign_from_swaps,
	factor_from_squares,
	geometric_prod_bits


@testset "bits <-> indices" begin
	for bits ∈ [0b000, 0b101, 0b110, 0b1101001, rand(UInt16, 3)..., rand(UInt32, 3)...]
		@test indices_to_bits(bits_to_indices(bits)) == bits
	end
end

@testset "bit permutations" begin
	@test first(BitPermutations{UInt8}(3)) === UInt8(0b111)
	@test first(BitPermutations{UInt32}(3)) === UInt32(0b111)
	@test collect(Iterators.take(BitPermutations(2), 3)) == [0b011, 0b101, 0b110]
	@test collect(Iterators.take(BitPermutations(5), 3)) == [0b011111, 0b101111, 0b110111]
end

@testset "bits <-> multivector index" begin
	cases = [
		(1, 0) => 0b000,
		(1, 1) => 0b001,
		(2, 1) => 0b010,
		(1, 2) => 0b011,
		(3, 1) => 0b100,
		(2, 2) => 0b101,
		(3, 2) => 0b110,
		(1, 3) => 0b111,
		(1, 8) => 0b01111_1111,
		(2, 8) => 0b10111_1111,
	]
	@testset "no. $(i) of grade $k: $bits" for ((i, k), bits) ∈ cases
		@test mv_index_to_bits(i, k) == bits
		@test bits_to_mv_index(bits) == i
	end

	for _ ∈ 1:10, bits ∈ (rand(UInt8), rand(UInt16))
		# warning: bits_to_mv_index() is typically slow for UInt32 or higher
		@test mv_index_to_bits(bits_to_mv_index(bits), count_ones(bits)) == bits
	end
end

@testset "bits <-> mixed multivector index" begin
	for dim ∈ 1:2:8, i ∈ 1:2:2^dim
		@test bits_to_mmv_index(mmv_index_to_bits(i, dim), dim) == i
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