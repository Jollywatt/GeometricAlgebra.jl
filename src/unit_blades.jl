grade(bits::Unsigned) = count_ones(bits)

bits_scalar() = unsigned(0)
bits_first_of_grade(k) = (unsigned(1) << k) - unsigned(1)

bits_basis_vector(i) = unsigned(1) << (i - 1)
# bits_has_index(bits, i) = isone(bits >> (i - 1) & 1)

function bits_to_indices(bits)
	indices = Int[]
	i = 1
	while bits > 0
		if isone(bits & 1)
			push!(indices, i)
		end
		i += 1
		bits >>= 1
	end
	indices
end

function indices_to_bits(indices)
	bits = bits_scalar()
	for i ∈ indices
		bits += bits_basis_vector(i)
	end
	bits
end


"""
	bits_to_indices(bits)

Return the positions of the ones in the unsigned integer `bits`.

Used to convert between representations of a unit blade.
See also [`$(repr(indices_to_bits))`](@ref).

Examples
===
```jldoctest
julia> GeometricAlgebra.bits_to_indices(0b1001101)
4-element Vector{Int64}:
 1
 3
 4
 7
```
"""
bits_to_indices

"""
	indices_to_bits(indices)

Create unsigned integer with bits at the positions given in the vector `indices`.

Used to convert between representations of a unit blade.
See also [`$(repr(bits_to_indices))`](@ref).

Examples
===
```jldoctest
julia> GeometricAlgebra.indices_to_bits([1, 2, 5]) |> UInt16 |> bitstring
"0000000000010011"
```

"""
indices_to_bits

# # this version is 10% slower
# function bits_to_linear_index(bits::Unsigned)
# 	i = 1
# 	for b ∈ FixedGradeBits(grade(bits))
# 		if b >= bits
# 			return i
# 		end
# 		i += 1
# 	end
# end
function bits_to_linear_index(bits::Unsigned)
	i = 1
	u = bits_first_of_grade(grade(bits))
	while u < bits
		u = next_bit_permutation(u)
		i += 1
	end
	i
end
function linear_index_to_bits(i, k)
	bits = bits_scalar()
	for b ∈ Iterators.take(FixedGradeBits(k), i)
		bits = b
	end
	bits
end

"""
Return the smallest uint larger than the one given which has
the same number of binary ones.
Algorithm is [Gosper's hack](http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation).

```
julia> GeometricAlgebra.next_bit_permutation(0b1011) |> bitstring
"00001101"
```
"""
function next_bit_permutation(v::Unsigned)
	I = one(v)
	t = v | (v - I)
	(t + I) | (((~t & -~t) - I) >> (trailing_zeros(v) + I))
end



struct FixedGradeBits
	grade::Int
end

function Base.iterate(fgb::FixedGradeBits, bits=bits_first_of_grade(fgb.grade))
	(bits, next_bit_permutation(bits))
end

Base.IteratorSize(::FixedGradeBits) = Base.IsInfinite()



"""
Compute sign flips of blade product due to transposing basis vectors.
(The full sign of the product will also depend on the basis norms.)
"""
function sign_from_swaps(a::Unsigned, b::Unsigned)
	swaps = 0
	while b > 0
		a >>= 1
		swaps += count_ones(a)*(b & 1)
		b >>= 1
	end
	iseven(swaps) ? 1 : -1
end


"""
Compute the overall factor arising from the geometric product between
repeated basis vectors.
"""
function factor_from_squares(sig, squares)
	factor = 1
	bv = 1
	while squares > 0
		if isone(squares & 1)
			factor *= sig[bv]
		end
		squares >>= 1
		bv += 1
	end
	factor
end


function geometric_prod_bits(sig, a::Unsigned, b::Unsigned)
	factor = sign_from_swaps(a, b)*factor_from_squares(sig, a & b)
	bits = a ⊻ b
	factor, bits
end