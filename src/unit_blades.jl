"""
Operations on bits representing `unit blades'. E.g., the
2-blade of unit norm `e₂∧e₃ ≡ e₂₃` is represented by `0b110`.
"""

grade(bits::Unsigned) = count_ones(bits)

bits_scalar() = unsigned(0)
bits_first_of_grade(k) = (unsigned(1) << k) - unsigned(1)

bits_basis_vector(i) = unsigned(1) << (i - 1)
# bits_has_index(bits, i) = isone(bits >> (i - 1) & 1)


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

struct FixedGradeBits
	grade::Int
end

function Base.iterate(fgb::FixedGradeBits, bits=bits_first_of_grade(fgb.grade))
	(bits, next_bit_permutation(bits))
end

Base.IteratorSize(::FixedGradeBits) = Base.IsInfinite()

"""
	bits_of_grade(k[, n])

Generate basis blade bits of grade `k` in ascending order.
Yields all basis blades in the dimension `n`, if given, otherwise iterate indefinitely.

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
bits_of_grade(k) = FixedGradeBits(k)
bits_of_grade(k, n) = Iterators.take(FixedGradeBits(k), binomial(n, k))



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

function mv_index_to_bits(ith, k)
	bits = bits_scalar()
	for b ∈ Iterators.take(bits_of_grade(k), ith)
		bits = b
	end
	bits
end


# const BINOMIAL_SUMS = Dict{Int,Vector{Int}}()
# function binomial_sum(n, k)
# 	k <= 0 && return 0
# 	n <= k && return 1<<n
# 	if !(n in keys(BINOMIAL_SUMS))
# 		BINOMIAL_SUMS[n] = cumsum(binomial.(n, 0:n - 1))
# 	end
# 	# @show BINOMIAL_SUMS
# 	BINOMIAL_SUMS[n][k]
# end

function multivector_index_offset(k, dim)
	ith = 0
	for i in 0:k - 1
		ith += binomial(dim, i)
	end
	ith
end
function bits_to_mmv_index(bits, dim)
	# memoized version only faster for rather large k (~10)
	# return 1 + binomial_sum(dim, grade(bits)) + bits_to_mv_index(bits)

	multivector_index_offset(grade(bits), dim) + bits_to_mv_index(bits)
end


const MULTIVECTOR_INDICES = Dict{Int,Vector{UInt}}()
function mmv_index_to_bits(ith, dim)
	if !(dim in keys(MULTIVECTOR_INDICES))
		bits = unsigned.(0:2^dim - 1)
		MULTIVECTOR_INDICES[dim] = bits[sortperm(bits_to_mmv_index.(bits, dim))]
	end
	MULTIVECTOR_INDICES[dim][ith]
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
Compute the overall factor arising from the geometric product between
repeated basis vectors.
"""
function factor_from_squares(sig, squares::Unsigned)
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

# geometric_square_sign(sig, a::Unsigned) = sign_from_swaps(a, a)*factor_from_squares(sig, a)
