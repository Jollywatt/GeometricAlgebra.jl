using GeometricAlgebra:
	bits_to_indices, indices_to_bits

@testset "bits <-> indices" begin
	for bits ∈ [0b000, 0b101, 0b110, 0b1101001, rand(UInt16), rand(UInt64)]
		# expect to fail for UInt128
		@test indices_to_bits(bits_to_indices(bits)) == bits
	end
end