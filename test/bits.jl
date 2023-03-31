using GeometricAlgebra:
	bits_to_indices,
	indices_to_bits,
	bits_of_grade,
	componentbits,
	componentindex,
	sign_from_swaps,
	factor_from_squares,
	geometric_prod_factor


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

@testset "componentbits" begin
	for n ∈ 0:8, k ∈ 0:n
		bits = collect(componentbits(n, k))
		@test length(bits) == binomial(n, k)
		@test all(count_ones.(bits) .== k)

		@test length(collect(componentbits(n, 0:n))) == 2^n

		if n > 1
			@test length(collect(componentbits(n, 0:2:n))) == 2^(n - 1)
			@test length(collect(componentbits(n, (0, n)))) == 2
		end
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
	@test factor_from_squares((0, 1, 1), 0b101) == 0
end

@testset "geometric_prod_factor()" begin
	@test geometric_prod_factor((1, 1, 1), 0b001, 0b010) == +1
	@test geometric_prod_factor((1, 1, 1), 0b010, 0b001) == -1
	@test geometric_prod_factor((1, 1, 1), 0b011, 0b001) == -1
	@test geometric_prod_factor((0, 1, 1), 0b001, 0b101) == 0
end