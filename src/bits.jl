"""
Operations on bits representing “unit blades”. E.g., the
2-blade of unit norm `e₂∧e₃ ≡ e₂₃` is represented by `0b110`.
"""

#= Conversion from Indices to Bits =#


"""
	bits_to_indices(bits)

Return the positions of the ones in the unsigned integer `bits`.

Used to convert between representations of a unit blade.
Inverse of [`indices_to_bits`](@ref).

# Examples
```jldoctest
julia> GeometricAlgebra.bits_to_indices(0b1001101)
4-element Vector{Int64}:
 1
 3
 4
 7
```
"""
function bits_to_indices(bits::Unsigned)
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
	indices_to_bits(indices, T=UInt)

Create unsigned integer with bits at the positions given in the vector `indices`.

Used to convert between representations of a unit blade.
Inverse of [`bits_to_indices`](@ref).

!!! warning
	Only correct if `maximum(indices)` does not exceed number of bits in `T`.

# Examples
```jldoctest
julia> GeometricAlgebra.indices_to_bits([1, 2, 5], UInt8) |> bitstring
"00010011"
```
"""
function indices_to_bits(indices, T=UInt)
	bits = zero(T)
	for i ∈ indices
		bits += one(T) << (i - 1)
	end
	bits
end


"""
Return the smallest uint larger than the one given which has
the same number of binary ones.
Algorithm is [Gosper’s hack](http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation).

# Examples
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
	n::Int
end
BitPermutations(n) = BitPermutations{UInt}(n)

function Base.iterate(itr::BitPermutations{T}, bits=(one(T) << itr.n) - one(T)) where T
	(bits, next_bit_permutation(bits))
end

Base.IteratorSize(::BitPermutations) = Base.IsInfinite()
Base.eltype(::Type{BitPermutations{T}}) where {T} = T

"""
	bits_of_grade(k[, dim])

Generate basis blade bits of grade `k` in ascending order.
Yields all basis blades in the dimension `dim`, if given, otherwise iterates indefinitely.

# Examples
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
bits_of_grade(k, n) = Iterators.take(BitPermutations(k), binomial(n, k))

bits_dual(n, bits) = first(bits_of_grade(n)) ⊻ bits



"""
	componentbits(n, k)
	componentbits(::Val{N}, ::Val{K})

Vector of bits corresponding to components of an `n`-dimensional
`Multivector` of grade(s) `k`.

Bits are ordered first by grade (`count_ones`),
then lexicographically (in ascending numerical order).

Passing `Val` arguments calls a faster, memoized method.

# Examples
```jldoctest; setup = :(using GeometricAlgebra: componentbits)
julia> componentbits(4, 2) .|> UInt8 .|> bitstring
6-element Vector{String}:
 "00000011"
 "00000101"
 "00000110"
 "00001001"
 "00001010"
 "00001100"

julia> componentbits(3, 0:3) .|> UInt8 .|> bitstring
8-element Vector{String}:
 "00000000"
 "00000001"
 "00000010"
 "00000100"
 "00000011"
 "00000101"
 "00000110"
 "00000111"
```
"""
componentbits(n, k) = collect(Iterators.flatten(bits_of_grade.(k, n)))
@generated componentbits(::Val{N}, ::Val{K}) where {N,K} = componentbits(N, K)




#= Geometric Products of Bits =#

"""
	sign_from_swaps(a::Unsigned, b::Unsigned)

Compute sign flips of blade product due to transposing basis vectors into sorted order.
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
repeated basis vectors.
"""
function factor_from_squares(sig, bits::Unsigned)
	factor = 1
	bv = 1
	while bits > 0
		if isone(bits & 1)
			factor *= basis_vector_norm(sig, bv)
		end
		bits >>= 1
		bv += 1
	end
	factor
end

geometric_prod_matrix(sig) = let bits = 0:bits_first_of_grade(dimension(sig))
	Int[sign_from_swaps(a, b)factor_from_squares(sig, a & b) for a ∈ bits, b ∈ bits]
end
@generated geometric_prod_matrix(::Val{Sig}) where {Sig} = geometric_prod_matrix(Sig)

@inline function geometric_prod_factor(sig, a, b)::Int
	if false # dimension(sig) <= 10 # limit geometric_prod_matrix cache to 8MB
		# doesn’t seem worth it...
		# TODO: proper benchmark test
		geometric_prod_matrix(Val(sig))[begin + a, begin + b]
	else
		sign_from_swaps(a, b)factor_from_squares(sig, a & b)
	end
end

"""
	geometric_prod_bits(sig, a::Unsigned, b::Unsigned)

Compute the geometric product between unit blades. Returns a tuple
of the overall scalar factor and the resulting unit blade.
"""
geometric_prod_bits(sig, a::Unsigned, b::Unsigned) = (geometric_prod_factor(sig, a, b), a ⊻ b)

"""
	reversion_sign(k) = mod(k, 4) <= 1 ? +1 : -1

Sign from reversing a ``k``-vector.
"""
reversion_sign(k) = mod(k, 4) <= 1 ? +1 : -1

geometric_square_factor(sig, a::Unsigned) = reversion_sign(count_ones(a))factor_from_squares(sig, a)
