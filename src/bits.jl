"""
Operations on bits representing “unit blades”. E.g., the
2-blade of unit norm `e₂∧e₃ ≡ e₂₃` is represented by `0b110`.
"""

bits_scalar() = unsigned(0)
bits_first_of_grade(T::Type{<:Unsigned}, k) = (one(T) << k) - one(T)
bits_first_of_grade(k) = bits_first_of_grade(UInt, k)
bits_basis_vector(i) = unsigned(1) << (i - 1)


"""
	bits_to_indices(bits)

Return the positions of the ones in the unsigned integer `bits`.

Used to convert between representations of a unit blade.
Inverse of [`indices_to_bits`](@ref).

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

"""
	indices_to_bits(indices)

Create unsigned integer with bits at the positions given in the vector `indices`.

Used to convert between representations of a unit blade.
Inverse of [`bits_to_indices`](@ref).

Examples
===
```jldoctest
julia> GeometricAlgebra.indices_to_bits([1, 2, 5]) |> UInt16 |> bitstring
"0000000000010011"
```

"""
function indices_to_bits(indices)
	bits = bits_scalar()
	for i ∈ indices
		bits += bits_basis_vector(i)
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



"""
	BitPermutations{T}(n)

Infinite iterator returning all unsigned integers of type `T`,
in ascending order, for which `Base.count_ones` is `n`.
"""
struct BitPermutations{T<:Unsigned}
	n::T
end
BitPermutations(n) = BitPermutations{UInt}(n)

function Base.iterate(fgb::BitPermutations{T}, bits=bits_first_of_grade(T, fgb.n)) where T
	(bits, next_bit_permutation(bits))
end

Base.IteratorSize(::BitPermutations) = Base.IsInfinite()



"""
	bits_of_grade(k[, dim])

Generate basis blade bits of grade `k` in ascending order.
Yields all basis blades in the dimension `dim`, if given, otherwise iterate indefinitely.

Example
===
```jldoctest
julia> GeometricAlgebra.bits_of_grade(2, 4) .|> UInt8 .|> bitstring
6-element Vector{String}:
 "00000011"
 "00000101"
 "00000110"
 "00001001"
 "00001010"
 "00001100"
```
"""
bits_of_grade(k) = BitPermutations(k)
bits_of_grade(k, dim) = Iterators.take(BitPermutations(k), binomial(dim, k))

"""
	bits_to_mv_index(bits::Unsigned)

Convert a unit blade `bits` to a linear index for accessing components of a `Multivector`. 
"""
function bits_to_mv_index(bits::Unsigned)
	# From combinatorial number systems, an explicit formula is:
	# sum(binomial(i, c[i]) for i in 1:length(c)) where c = bits_to_indices(bits)
	ith = 1
	c = 0
	k = 1
	while bits > 0
		if isone(bits & 1)
			ith += binomial(c, k)
			k += 1
		end
		bits >>= 1
		c += 1
	end
	ith
end

"""
	mv_index_to_bits(ith, k)

Convert a linear index for grade-`k` `Multivector` components
into the corresponding unit blade.
"""
function mv_index_to_bits(ith, k)
	bits = bits_scalar()
	for b ∈ Iterators.take(bits_of_grade(k), ith)
		bits = b
	end
	bits
end


function multivector_index_offset(k, dim)
	ith = 0
	for i in 0:k - 1
		ith += binomial(dim, i)
	end
	ith
end

"""
	bits_to_mmv_index(bits::Unsigned)

Convert a unit blade `bits` to a linear index for accessing components of a `MixedMultivector`. 
"""
function bits_to_mmv_index(bits, dim)
	multivector_index_offset(count_ones(bits), dim) + bits_to_mv_index(bits)
end

# range of MixedMultivector corresponding to the grade k components
mmv_slice(k, dim) = multivector_index_offset(k, dim) .+ (1:binomial(dim, k))


const MULTIVECTOR_INDICES = Dict{Int,Vector{UInt}}()

"""
	mmv_index_to_bits(ith, dim)

Convert a linear index for `MixedMultivector` components into the corresponding
unit blade. Uses a lookup table `$(@__MODULE__).MULTIVECTOR_INDICES` for speed.
"""
function mmv_index_to_bits(ith, dim)
	mmv_index_to_bits(dim)[ith]
end
function mmv_index_to_bits(dim)
	if !(dim in keys(MULTIVECTOR_INDICES))
		bits = unsigned.(0:2^dim - 1)
		MULTIVECTOR_INDICES[dim] = bits[sortperm(bits_to_mmv_index.(bits, dim))]
	end
	MULTIVECTOR_INDICES[dim]
end

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
	factor_from_squares(sig, bits::Unsigned)

Compute the overall factor arising from the geometric product between
repeated basis vectors
"""
function factor_from_squares(sig, bits::Unsigned)
	factor = 1
	bv = 1
	while bits > 0
		if isone(bits & 1)
			factor *= sig[bv]
		end
		bits >>= 1
		bv += 1
	end
	factor
end


"""
	geometric_prod_bits(sig, a::Unsigned, b::Unsigned)

Compute the geometric product between unit blades. Returns a tuple
of the overall scalar factor and the resulting unit blade.
"""
function geometric_prod_bits(sig, a::Unsigned, b::Unsigned)
	factor = sign_from_swaps(a, b)*factor_from_squares(sig, a & b)
	bits = a ⊻ b
	factor, bits
end

reversion_sign(k) = mod(k, 4) <= 1 ? +1 : -1
geometric_square_factor(sig, a::Unsigned) = reversion_sign(count_ones(a))*factor_from_squares(sig, a)