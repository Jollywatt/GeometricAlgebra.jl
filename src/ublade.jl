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
sequence of basic vectors, only the last representation is possible.

One this project matures enough to think about performance, it could be
decided whether BitVectors offer any advantage.
Practically, I can only see uints being used, and otherwise vectors of symbols.

=#

# for debugging:
Base.show(io::IO, ::MIME"text/plain", a::Unsigned) = print(io, join(reverse(digits(a, base=2; pad=8sizeof(a)))))

# first grade-k unit blade basis under lexicographic ordering; e.g., 0b0111, [1,2,3]
ublade_first_of_grade(T::Type{<:Unsigned}, k) = one(T) << k - one(T)
ublade_first_of_grade(T::Type{<:Vector}, k) = T(1:k)
ublade_first_of_grade(T::Type{<:BitVector}, k) = trues(k)

# scalar unit blade basis of the given type; e.g., 0b0000, []
ublade_scalar(T, sig) = ublade_first_of_grade(T, 0)

# pseudoscalar unit blade basis of given type
ublade_vol(T, sig) = ublade_first_of_grade(T, dim(sig))

# unit blade basis of ith basis vector 
ublade_bv(T::Type{<:Unsigned}, i) = one(T) << (i - 1)
ublade_bv(T::Type{Vector{<:Integer}}, i) = [i]
ublade_bv(T::Type{BitVector}, i) = [falses(i - 1); true]

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
	for i ∈ eachindex(uv)
		uv[i] = get(u, i, 0) ⊻ get(v, i, 0)
	end
	uv
end



# UNIT BLADE PROMOTION & CONVERSION

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


# this might be simplified if BitVectors are banished
promote_ublade_type(as::Type{<:Unsigned}...) = promote_type(as...)
for T ∈ [Vector{<:Integer}, BitVector]
	@eval promote_ublade_type(a::Type{<:Unsigned}, b::Type{<:$T}) = b
	@eval promote_ublade_type(a::Type{<:$T}, b::Type{<:Unsigned}) = a
	@eval promote_ublade_type(a::Type{<:$T}, b::Type{<:$T}) = promote_type(a, b)
end

promote_ublade_type(::Type{T}...) where T = T

function promote_ublades(ublades...)
	T = promote_ublade_type(typeof.(ublades)...)
	tuple(convert_ublade.(T, ublades)...)
end





# UNIT BLADE OPERATIONS


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
Give the next unit blade of the same grade in lexicographic order.
"""
next_ublade(u) = next_ublade!(copy(u))

next_ublade!(u::Unsigned) = next_bit_permutation(u)
function next_ublade!(u::Vector{<:Integer})
	for i ∈ firstindex(u):lastindex(u) - 1
		if u[i + 1] <= u[i] + 1
			u[i] = i
			i += 1
		else
			u[i] = u[i] + 1
			return u
		end
	end
	u[end] += 1
	u
end
function next_ublade!(u::BitVector)
	nwrapped = 1
	for i ∈ 1:length(u) - 1
		u[i] || continue
		u[i] = false
		if u[i + 1]
			u[nwrapped] = true
			nwrapped += 1
		else
			u[i + 1] = true
			return u
		end
	end
	u[end] = false
	push!(u, true)
	u
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
function lindex2ublade(T, k, i)
	ublade = ublade_first_of_grade(T, k)
	for _ ∈ 2:i
		ublade = next_ublade!(ublade)
	end
	ublade
end

function ublade2lindex(ublade::Unsigned)
	#TODO: does this have a better algorithm?
	lindex = 1
	u = ublade_first_of_grade(typeof(ublade), ublade_grade(ublade))
	while u < ublade
		u = next_ublade(u)
		lindex += 1
	end
	lindex
end
function ublade2lindex(ublade::BitVector)
	lindex = 1
	k = 0
	for (cₖ, bit) ∈ enumerate(ublade)
		if bit
			k += 1
			lindex += binomial(cₖ - 1, k)
		end
	end
	lindex
end
function ublade2lindex(ublade::Vector{<:Integer})
	lindex = 1
	for (k, cₖ) ∈ enumerate(ublade)
		lindex += binomial(cₖ - 1, k)
	end
	lindex
end



"""
	FixedGradeBlades{B}(k, n)

Generates all grade `k` unit blades in `n` dimensions ordered lexicographically.
Unit blades of type `B` are produced.
"""
struct FixedGradeBlades{B}
	grade::Int
	dimension::Int
end
function Base.iterate(fgb::FixedGradeBlades{B}) where B
	remaining = binomial(fgb.dimension, fgb.grade)
	if remaining > 0
		ublade = ublade_first_of_grade(B, fgb.grade)
		(ublade, (ublade, remaining - 1))
	else
		nothing
	end
end
function Base.iterate(fgb::FixedGradeBlades{B}, (ublade, remaining)) where B
	if remaining > 0
		ublade = next_ublade(ublade)
		(ublade, (ublade, remaining - 1))
	else
		nothing
	end
end
Base.length(fgb::FixedGradeBlades) = binomial(fgb.dimension, fgb.grade)
Base.eltype(::FixedGradeBlades{B}) where B = B



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
function ubladeprod(sig, bvs::Vector{<:Integer}; startfrom=0)
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
function _ubladeprod(sig, a, b)
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
ubladeprod(sig, a::Unsigned, b::Unsigned) = _ubladeprod(sig, a, b)
ubladeprod(sig, a::BitVector, b::BitVector) = _ubladeprod(sig, a, b)
ubladeprod(sig, a::Union{<:Unsigned,BitVector}) = (1, a)


