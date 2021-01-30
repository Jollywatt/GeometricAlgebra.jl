#= UNIT BLADES

Unit blades represent the wedge product of orthonormal
basis vectors such as ``v_1∧v_3∧v_4`` and may be represented as:
 - Unsigned           0b1101
 - Vector{<:Integer}  [1, 3, 4]
 - Vector{<:Symbol}   [:v1, :v3, :v4]
depending on the multivector type and metric signature.
This module tries to stay agnostic with respect to ublade type.
However, in signatures of unspecified dimension which lack a canonical
sequence of basic vectors, only the last representation is possible.

Homogeneous `k`-multivectors have a canonical lexicographical ordering of basis blades.
E.g., the 6 basis 2-blades in 4 dimensions are, in order:
order	1		2		3		4		5		6
uints	0011	0101	0110	1001	1010	1100
vects	[1, 2]	[1, 3]	[2, 3]	[1, 4]	[2, 4]	[3, 4]
Viewed as int vects, k-blades form a combinatorial number system,
where their value is the 'order' index and their expression in
the combinatorial number system is the corresponding `vect`.
https://en.wikipedia.org/wiki/Combinatorial_number_system
=#


# UNIT BLADE UTILITY FUNCTIONS

# first grade-k unit blade basis under lexicographic ordering; e.g., 0b0111, [1,2,3]
ublade_first_of_grade(T::Type{<:Unsigned}, k) = one(T) << k - one(T)
ublade_first_of_grade(T::Type{<:Vector{<:Integer}}, k) = T(1:k)

# scalar unit blade basis of the given type; e.g., 0b0000, []
ublade_scalar(T) = ublade_first_of_grade(T, 0)
ublade_scalar(::Type{<:Vector{Symbol}}) = Symbol[]

ublade_vol(sig, T) = ublade_first_of_grade(T, dim(sig))

# unit blade basis of ith basis vector 
ublade_bv(T::Type{<:Unsigned}, i) = one(T) << (i - 1)
ublade_bv(T::Type{Vector{<:Integer}}, i) = [i]

ublade_has_bv(ublade::Unsigned, i) = !iszero(ublade & 1 << (i - 1))
ublade_has_bv(ublade::Vector, i) = i ∈ ublade

ublade_grade(ublade::Unsigned) = count_ones(ublade)
ublade_grade(ublade::Vector) = length(ublade)

ublade_xor(u::Unsigned, v::Unsigned) = u ⊻ v
ublade_xor(u::Vector, v::Vector) = symdiff(u, v)

# iterate indices of basis vectors present in unit blade
ublade_bvs(sig, ublade) = convert_ublade(sig, Vector{Int}, ublade)

# scalar square of unit blade
ublade_square(sig, ublade) = prod(sig[i] for i ∈ ublade_bvs(sig, ublade))



# UNIT BLADE PROMOTION & CONVERSION
# convert_ublade assumes ublades in canonical form
# e.g., both [1, 2] and [2, 1] are converted to 0b011

convert_ublade(sig, ::Type{T}, ublade::T) where T = ublade

# Unsigned <-> Vector{<:Integer}
function convert_ublade(sig, T::Type{<:Unsigned}, ublade::Vector{<:Integer})
	@assert issorted(ublade)
	u = zero(T)
	for bv ∈ ublade
		u += one(T) << (bv - 1)
	end
	u
end
function convert_ublade(sig, T::Type{<:Vector{<:Integer}}, ublade::Unsigned)
	T(filter(i -> ublade_has_bv(ublade, i), 1:8sizeof(ublade)))
end

convert_ublade(sig, T::Type{<:Unsigned}, ublade::Unsigned) = convert(T, ublade)

# Vector{<:Integer} <-> Vector{Symbol}
convert_ublade(sig, T::Type{<:Vector{<:Integer}}, ublade::Vector{Symbol}) =
	T([let i = findfirst(==(symbol), signature_labels(sig))
		isnothing(i) && error("signature $sig has no label $(repr(symbol))")
		i
	end for symbol ∈ ublade])
convert_ublade(sig, T::Type{<:Vector{Symbol}}, ublade::Vector{<:Integer}) =
	T([signature_labels(sig)[i] for i ∈ ublade])

# Unsigned <-> Vector{Symbol}
convert_ublade(sig, T::Type{<:Unsigned}, ublade::Vector{Symbol}) =
	convert_ublade(sig, T, convert_ublade(sig, Vector{signed(T)}, ublade))
convert_ublade(sig, T::Type{<:Vector{Symbol}}, ublade::Unsigned) =
	convert_ublade(sig, T, convert_ublade(sig, Vector{signed(typeof(ublade))}, ublade))

"""
	convert_ublade(a::AbstractMultivector, ublade)

Convert unit blade `ublade` to the representation employed by `a`.
"""

# first argument should be a::AbstractMultivector
convert_ublade(a, ublade) = convert_ublade(signature(a), keytype(a), ublade)

# fallback
# convert_ublade(sig, ::Type{T}, ublade::S) where {T,S} = convert(T, ublade)


promote_ublade_type(as::Type{<:Unsigned}...) = promote_type(as...)
promote_ublade_type(a::Type{<:Unsigned}, b::Type{<:Vector{<:Integer}}) = b
promote_ublade_type(a::Type{<:Vector{<:Integer}}, b::Type{<:Unsigned}) = a
promote_ublade_type(a::Type{<:Vector{<:Integer}}, b::Type{<:Vector{<:Integer}}) = promote_type(a, b)
promote_ublade_type(::Type{T}...) where T = T

function promote_ublades(sig, ublades...) # used only by ==(::Blade, ::Blade) so far
	T = promote_ublade_type(typeof.(ublades)...)
	Tuple(convert_ublade(sig, T, u) for u ∈ ublades)
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


"""
Convert linear index `i` of into unit `k`-blade.
The linear index (lindex) gives the index of the unit blade
sorted lexicographically.

```
julia> lindex2ublade.(UInt8, 3, 1:5) .|> bitstring
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
	lindex = 1
	u = ublade_first_of_grade(typeof(ublade), ublade_grade(ublade))
	while u < ublade
		u = next_ublade(u)
		lindex += 1
	end
	lindex
end
function ublade2lindex(ublade::Vector{<:Integer})
	lindex = 1
	for (k, c_k) ∈ enumerate(ublade)
		lindex += binomial(c_k - 1, k)
	end
	lindex
end



"""
	FixedGradeBlades{B}(k, n)

Generates all grade `k` unit blades in `n` dimensions ordered lexicographically.
Unit blades of type `B` are produced.

Examples
===
```jldoctest; setup = :( using GeometricAlgebra: FixedGradeBlades )
julia> FixedGradeBlades{UInt8}(2, 3) .|> bitstring
3-element Array{String,1}:
 "00000011"
 "00000101"
 "00000110"

julia> FixedGradeBlades{Vector{Int}}(3, 4) |> collect
4-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]
```
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
function ubladeprod(sig, bvs::Vector; startfrom=0)
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
ubladeprod(sig, lbases, rbases) = ubladeprod(sig, vcat(lbases, rbases); startfrom=length(lbases))

"""
Compute sign flips of blade product due to transposing basis vectors.
(The full sign of the product will also depend on the basis norms.)
"""
function bladeprodsign(a::Unsigned, b::Unsigned)
	swaps = 0
	while b > 0
		a >>= 1
		swaps += count_ones(a)*(b & 1)
		b >>= 1
	end
	iseven(swaps)
end

# faster implementation
function ubladeprod(sig, a::Unsigned, b::Unsigned)
	factor = bladeprodsign(a, b) ? 1 : -1
	ab = a ⊻ b # generic on unsigned ints / bit vectors
	squares = a & b
	bv = 1
	while squares > 0
		if !iszero(squares & 1)
			factor *= sig[bv]
		end
		squares >>= 1
		bv += 1
	end
	factor, ab
end
ubladeprod(sig, a::Unsigned) = (1, a)

