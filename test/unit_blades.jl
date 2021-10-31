using GeometricAlgebra:
	bits_to_indices, indices_to_bits,
	bits_to_mv_index, mv_index_to_bits,
	bits_to_mmv_index, mmv_index_to_bits,
	sign_from_swaps, factor_from_squares

@testset "bits <-> indices" begin
	for bits ∈ [0b000, 0b101, 0b110, 0b1101001, rand(UInt16), rand(UInt64)]
		# expect to fail for UInt128
		@test indices_to_bits(bits_to_indices(bits)) == bits
	end
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
		@test mv_index_to_bits(bits_to_mv_index(bits), grade(bits)) == bits
	end
end

@testset "bits <-> mixed multivector index" begin
	for dim ∈ 1:8, i ∈ 1:2^dim
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
	@test factor_from_squares((1, 1, 1), 0b000) == 1
	@test factor_from_squares((1, 1, 1), 0b111) == 1
	@test factor_from_squares((-1, 1, 1), 0b001) == -1
	@test factor_from_squares((-1, 1, 1), 0b010) == +1
	@test factor_from_squares((-1, 1, 1), 0b111) == -1
	@test factor_from_squares((-1, -1, -1), 0b101) == +1
	@test factor_from_squares((-1, -1, -1), 0b111) == -1
end