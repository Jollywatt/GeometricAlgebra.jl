#= MULTIVECTOR TYPES

Type hierarchy:

           AbstractMultivector
            /               \
HomogeneousMultivector     MixedMultivector
   /        \                        
Blade   Multivector                 │
                                    │
           ╰─ CompositeMultivector ─╯

- `Blade`: a scalar multiple of a wedge product of orthogonal basis vectors.
   Note that the mathematical definition of a ``k``-blade is the wedge product
   of ``k`` different _vectors_, not necessarily basis vectors. Thus, not all
   ``k``-blades are representable as a `Blade` (but always as a `Multivector`).
- `Multivector`: a homogeneous multivector; a sum of same-grade blades.
- `MixedMultivector`: an inhomogeneous multivector. All elements in a geometric
   algebra can be represented as this type.

=#

"""
	AbstractMultivector{sig}

Supertype of all elements in the geometric algebra over the vector space with
metric signature `sig`, retrieved with the `signature` method.
"""
abstract type AbstractMultivector{sig} <: Number end


"""
	HomogeneousMultivector{sig,k} <: AbstractMultivector{sig}

Supertype of grade `k` elements in the geometric algebra with metric signature `sig`.
"""
abstract type HomogeneousMultivector{sig,k} <: AbstractMultivector{sig} end


"""
	Blade{sig,k,T} <: HomogeneousMultivector{sig,k}

A blade of grade `k` with basis blade `bits` and scalar coefficient of type `T`.

Parameters
===
- `sig`: metric signature defining the parent geometric algebra
- `k`: grade of the blade, equal to `grade(bits) === count_ones(bits)`
- `T`: type of the scalar coefficient
"""
struct Blade{sig,k,T} <: HomogeneousMultivector{sig,k}
	coeff::T
	bits::UInt
end
Blade{sig}(coeff::T, bits) where {sig,T} = Blade{sig,grade(bits),T}(coeff, bits)
Blade{sig,k}(coeff::T, bits) where {sig,k,T} = Blade{sig,k,T}(coeff, bits)
# warning: above method assumes `k == grade(bits)` without checking
bitsof(a::Blade) = a.bits

"""
	Multivector{sig,k,S} <: HomogeneousMultivector{sig,k}

A homogeneous multivector of grade `k` with storage type `S`.

Parameters
===
- `sig`: metric signature defining the parent geometric algebra
- `k`: grade of the multivector
- `S`: type in which the multivector components are stored; usually a vector-like or dictionary-like type
"""
struct Multivector{sig,k,S} <: HomogeneousMultivector{sig,k}
	components::S
end
Multivector{sig,k}(comps::S) where {sig,k,S} = Multivector{sig,k,S}(comps)


"""
	MixedMultivector{sig,S} <: AbstractMultivector{sig}

A (possibly inhomogeneous) multivector.

All elements of a geometric algebra are representable as a `MixedMultivector`.

Parameters
===
- `sig`: metric signature defining the parent geometric algebra
- `S`: type in which the multivector components are stored; usually a vector-like or dictionary-like type
"""
struct MixedMultivector{sig,S} <: AbstractMultivector{sig}
	components::S
end
MixedMultivector{sig}(comps::S) where {sig,S} = MixedMultivector{sig,S}(comps)


const CompositeMultivector{S} = Union{Multivector{sig,k,S},MixedMultivector{sig,S}} where {sig,k}
# const Scalar = Union{filter(T -> !(T <: AbstractMultivector), subtypes(Number))...}
const Scalar = Number
const OrType{X} = Union{X,Type{X}}



signature(::OrType{<:AbstractMultivector{sig}}) where sig = sig
dimension(a::OrType{<:AbstractMultivector}) = dimension(signature(a))
grade(::OrType{<:HomogeneousMultivector{sig,k} where sig}) where k = k


ncomponents(a::OrType{<:Multivector}) = binomial(dimension(a), grade(a))
ncomponents(a::OrType{<:MixedMultivector}) = 2^dimension(a)


Base.broadcastable(a::AbstractMultivector) = Ref(a)
Base.length(::AbstractMultivector) = error("$length is not defined for multivectors. Do you mean $(repr(ncomponents))()?")


"""
	eltype(a)

The numerical type of the components of a multivector or multivector type.
"""
Base.eltype(::Type{<:Blade{sig,k,T} where {sig,k}}) where T = T
Base.eltype(::Type{<:CompositeMultivector{S}}) where S = valtype(S)




#= ZERO & ONE CONSTRUCTORS =#

zeroslike(::Type{Vector{T}}, I...) where T = zeros(T, I...)
zeroslike(::Type{SparseVector{Tv,Ti}}, I...) where {Tv,Ti} = spzeros(Tv, Ti, I...)
zeroslike(::Type{S}, I...) where {S<:StaticVector} = zeros(S)

Base.zero(::Type{<:Blade{sig,k,T}}) where {sig,k,T} = Blade{sig,k,T}(zero(T), 0)
Base.zero(::Type{<:Multivector{sig,k,S}}) where {sig,k,S} = Multivector{sig,k}(zeroslike(S, binomial(dimension(sig), k)))
Base.zero(::Type{<:MixedMultivector{sig,S}}) where {sig,S} = MixedMultivector{sig}(zeroslike(S, 2^dimension(sig)))

# needed to override default implementation for `<:Number`
Base.zero(::T) where {T<:AbstractMultivector} = zero(T)

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector{<:AbstractVector}) = iszero(a.components)

Base.one(a::Type{<:Blade}) = Blade{signature(a)}(one(eltype(a)), bits_scalar())

Base.oneunit(a::Type{<:Blade}) = Blade{signature(a)}(one(eltype(a)), bitsof(a))
Base.oneunit(::T) where {T<:AbstractMultivector} = oneunit(T)




#= COMPONENT ACCESS =#

bits_to_index(a::OrType{<:Multivector}, bits) = bits_to_mv_index(bits)
bits_to_index(a::OrType{<:MixedMultivector}, bits) = bits_to_mmv_index(bits, dimension(a))
index_to_bits(a::OrType{<:Multivector}, ith) = mv_index_to_bits(ith, grade(a))
index_to_bits(a::OrType{<:MixedMultivector}, ith) = mmv_index_to_bits(ith, dimension(a))

_blades(a::Blade) = ((bitsof(a), a.coeff),)
_blades(a::Multivector) = zip(bits_of_grade(grade(a)), a.components)
_blades(a::MixedMultivector) = zip(mmv_index_to_bits.(1:ncomponents(a), dimension(a)), a.components)


blades(a::Blade) = (a,)
blades(a::Multivector{sig,k,<:AbstractVector}) where {sig,k} = (Blade{sig,k}(coeff, bits) for (coeff, bits) ∈ zip(a.components, bits_of_grade(k)))
blades(a::MixedMultivector{sig,<:AbstractVector}) where sig = (Blade{sig}(coeff, index_to_bits(a, i)) for (i, coeff) ∈ enumerate(a.components))


"""
	getcomponent(a::AbstractMultivector, I...)
	a[I...]

Return the component of the multivector `a` corresponding to the basis blade
defined by `I`.

Examples
===
```jldoctest
julia> v = basis((1, 1, 1));

julia> a = 1 + 2v[1] + 3v[2]v[3];

julia> a[], a[1], a[2], a[2,3]
(1, 2, 0, 3)

julia> GeometricAlgebra.getcomponent(a, 3, 2)
-3
```
"""
getcomponent
 
parity_sign(I) = iseven(parity(sortperm(collect(I)))) ? +1 : -1

function getcomponent(a::Blade, I...)
	indices_to_bits(I) == bitsof(a) || return zero(eltype(a))
	s = parity_sign(I)
	s*a.coeff
end

getcomponent(a::MixedMultivector) = first(a.components)

function getcomponent(a::MixedMultivector, I...)
	s = parity_sign(I)
	s*a.components[bits_to_index(a, indices_to_bits(I))]
end
function getcomponent(a::Multivector, I...)
	length(I) == grade(a) || return zero(eltype(a))
	s = parity_sign(I)
	s*a.components[bits_to_index(a, indices_to_bits(I))]
end


for T ∈ (Integer, Symbol)
	@eval Base.getindex(a::AbstractMultivector{sig}, I::$T...) where sig = getcomponent(a, normalize_bv_index.(Ref(sig), I)...)
end
Base.getindex(a::AbstractMultivector) = getcomponent(a) # for ambiguity

#= experimental notations for grade selection
What about a[(i,j)] for components and a[k] for grade projection?
Then if length(a) is the dimension, then a[end] is the pseudoscalar
Or a[grade=k], a[grade=end]
=#
# Base.getindex(a::AbstractMultivector, I::Vector{<:Integer}) = grade(a, Iterators.only(I))
# Base.getindex(a::AbstractMultivector; grade) = GeometricAlgebra.grade(a, grade)


