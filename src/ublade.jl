#= UNIT BLADES

Unit blades represent the wedge product of orthonormal
basis vectors such as v1∧v3∧v4 and may be represented as:
 - Unsigned           0b1101
 - BitVector          BitVector([1, 1, 0, 1])
 - Vector{<:Integer}  [1, 3, 4]
 - Vector{<:Symbol}   [:v1, :v3, :v4]
depending on the multivector type and metric signature.
This module tries to stay agnostic with respect to ublade type.
However, in signatures of unspecified dimension which lack a canonical
ordering of basic vectors, only the last representation is possible.

One this project matures enough to think about performance, it could be
decided whether BitVectors offer any advantage.
Practically, I can only see uints being used, and otherwise vectors of symbols.

=#

# scalar unit blade basis of the given type
ublade_scalar(sig, T::Type{<:Unsigned}) = zero(T)
ublade_scalar(sig, T::Type{<:Vector}) = T()
ublade_scalar(sig, ::Type{<:BitVector}) = falses(dim(sig))

# pseudoscalar unit blade basis of given type
ublade_vol(sig, T::Type{<:Unsigned}) = one(T) << dim(sig) - one(T)
ublade_vol(sig, T::Type{<:Vector}) = T(1:dim(sig))
ublade_vol(sig, ::Type{<:BitVector}) = trues(dim(sig))

# unit blade basis of ith basis vector 
ublade_bv(sig, T::Type{<:Unsigned}, i) = one(T) << (i - 1)
ublade_bv(sig, T::Type{Vector{<:Integer}}, i) = [i]
ublade_bv(sig, T::Type{BitVector}, i) = BitVector(1:dim(sig) .== i)

ublade_grade(ublade::Unsigned) = count_ones(ublade)
ublade_grade(ublade::Vector) = length(ublade)
ublade_grade(ublade::BitVector) = sum(ublade)

has_bv(ublade::Unsigned, i) = !iszero(ublade & 1 << (i - 1))
has_bv(ublade::Vector, i) = i ∈ ublade
has_bv(ublade::BitVector, i) = ublade[i]

ublade_xor(u::Unsigned, v::Unsigned) = u ⊻ v
ublade_xor(u::Vector, v::Vector) = symdiff(u, v)
function ublade_xor(u::BitVector, v::BitVector)
	uv = falses(max(length(u), length(v)))
	for i ∈ 1:length(uv)
		uv[i] = get(u, i, 0) ⊻ get(v, i, 0)
	end
	uv
end



# UNIT BLADE CONVERSIONS

function convert_ublade(T::Type{<:Unsigned}, ublade::Vector{<:Integer})
	u = zero(T)
	for bv ∈ ublade
		u += one(T) << (bv - 1)
	end
	u
end
function convert_ublade(T::Type{<:Vector{<:Integer}}, ublade::Unsigned)
	T(filter(i -> has_bv(ublade, i), 1:8sizeof(ublade)))
end

convert_ublade(T::Type{<:Unsigned}, ublade::Unsigned) = convert(T, ublade)
convert_ublade(::Type{T}, ublade::T) where T = ublade

# fallback
convert_ublade(::Type{T}, ublade::S) where {T,S} = convert(T, ublade)

#TODO: implement rest -- if they're really needed?



"""
Return the smallest uint larger than the one given which has
the same number of binary ones.
Algorithm is Gosper's hack, documented at
http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation.

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
Convert linear index `i` of unit `k`-blade into unsigned ublade form.
The `i`th ublade of grade `k` is the `i`th uint with `k` binary ones,
sorted lexicographically / in ascending order.

```
julia> GeometricAlgebra.lindex2ublade.(3, 1:5) .|> UInt8 .|> bitstring
5-element Vector{String}:
 "00000111"
 "00001011"
 "00001101"
 "00001110"
 "00010011"
```
"""
function lindex2ublade(k, i)
	I = unsigned(1)
	ublade = I << k - I
	while i > 1
		ublade = next_bit_permutation(ublade)
		i -= 1
	end
	ublade
end

function ublade2lindex(ublade)
	#TODO: does this have a better algorithm?
	i = 1
	u = lindex2ublade(ublade_grade(ublade), 1)
	while u < ublade
		u = next_bit_permutation(u)
		i += 1
	end
	i
end


# GEOMETRIC PRODUCT OF UNIT BLADES

"""
Compute the geometric product of unit basis vectors,
returning a scalar coefficient and the resulting unit blade.

```
julia> GeometricAlgebra.ubladeprod((1,1,1), [1, 2, 1, 2, 3])
(-1, [3])

julia> GeometricAlgebra.ubladeprod((1,1,1), 0b011, 0b111)
(-1, 0x04)
```
"""
function ubladeprod(sig, bvs; startfrom=0)
	sorted_bvs = bvs[1:startfrom]
	coeff = 1
	for b ∈ bvs[startfrom + 1:end]
		i = length(sorted_bvs)
		while i > 0
			if b < sorted_bvs[i] # basis vector belongs to the left
				coeff *= -1 # from anticommuting transposition
			elseif sorted_bvs[i] == b # repeated basis vector
				popat!(sorted_bvs, i) # basis vectors annihilate each other
				coeff *= sig[b] # and contribute square norm to overall coeff
				break
			else
				insert!(sorted_bvs, i + 1, b)
				break
			end
			i -= 1
		end
		if iszero(i)
			insert!(sorted_bvs, 1, b)
		end
	end
	coeff, sorted_bvs
end

# optimisation for partially-sorted prod
# only correct if `lbases` is in canonical form (strictly increasing)
ubladeprod(sig, lbases, rbases) = ubladeprod(sig, [lbases..., rbases...]; startfrom=length(lbases))

"""
Compute sign flips of blade product due to transposing basis vectors.
(The full sign of the product will also depend on the basis norms.)
"""
function bladeprodsign(a::BitVector, b::BitVector)
	iseven(sum((a << 1) .* cumsum(b))) ? 1 : -1
end
function bladeprodsign(a::Unsigned, b::Unsigned)
	swaps = 0
	while b > 0
		a >>= 1
		swaps += count_ones(a)*(b & 1)
		b >>= 1
	end
	iseven(swaps) ? 1 : -1
end

# faster implementation
function ubladeprod(sig, a::Unsigned, b::Unsigned)
	factor = bladeprodsign(a, b)
	ab = a .⊻ b # generic on unsigned ints / bit vectors
	squares = a .& b
	for i ∈ 1:dim(sig)
		if has_bv(squares, i)
			factor *= sig[i]
		end
	end
	factor, ab
end
ubladeprod(sig, a::Unsigned) = (1, a)


