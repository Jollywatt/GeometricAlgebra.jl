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

signature(::Type{<:AbstractMultivector{sig}}) where sig = sig
signature(::T) where {T<:AbstractMultivector} = signature(T)

dimension(::Union{T,Type{T}}) where {T<:AbstractMultivector} = dimension(signature(T))

Base.broadcastable(a::AbstractMultivector) = Ref(a)

"""
	HomogeneousMultivector{sig,k} <: AbstractMultivector{sig}

Supertype of grade `k` elements in the geometric algebra with metric signature `sig`.
"""
abstract type HomogeneousMultivector{sig,k} <: AbstractMultivector{sig} end

grade(::Type{<:HomogeneousMultivector{sig,k} where sig}) where k = k
grade(::T) where {T<:HomogeneousMultivector} = grade(T)


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

# for internal use only (assumes `k == grade(bits)` without checking)
Blade{sig,k}(coeff::T, bits) where {sig,k,T} = Blade{sig,k,T}(coeff, bits)
# (Blade{sig,k,T} where k)(coeff, bits) where {sig,T} = Blade{sig,grade(bits),T}(coeff, bits)
# Blade{sig,k,bits}(coeff::T) where {sig,k,bits,T} = Blade{sig,k,bits,T}(coeff)

"""
	bitsof(a::Blade{sig,k,bits,T}) = bits

The bits parameter of a blade, representing its basis `k`-blade.
"""
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

An generally inhomogeneous multivector.

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


# const UnionMultivectorTypes = Union{Blade,Multivector,AbstractMultivector}
const CompositeMultivector{S} = Union{Multivector{sig,k,S},MixedMultivector{sig,S}} where {sig,k}
# const Scalar = Union{filter(T -> !(T <: AbstractMultivector), subtypes(Number))...}
const Scalar = Number


"""
	eltype(a)

The numerical type of the components of a multivector or multivector type.
"""
Base.eltype(::Type{<:Blade{sig,k,T} where {sig,k}}) where T = T
Base.eltype(::Type{<:CompositeMultivector{S}}) where S = valtype(S)


ncomponents(::Type{<:Multivector{sig,k}}) where {sig,k} = binomial(dimension(sig), k)
ncomponents(::Type{<:MixedMultivector{sig}}) where sig = 2^dimension(sig)
ncomponents(::T) where {T<:AbstractMultivector} = ncomponents(T)
Base.length(::AbstractMultivector) = error("$length is not defined for multivectors. Do you mean $(repr(ncomponents))()?")


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





bits_to_index(a::Type{<:Multivector}, bits) = bits_to_mv_index(bits)
bits_to_index(a::Type{<:MixedMultivector}, bits) = bits_to_mmv_index(bits, dimension(a))
index_to_bits(a::Type{<:Multivector}, ith) = mv_index_to_bits(ith, grade(a))
index_to_bits(a::Type{<:MixedMultivector}, ith) = mmv_index_to_bits(ith, dimension(a))

bits_to_index(a::AbstractMultivector, b) = bits_to_index(typeof(a), b)
index_to_bits(a::AbstractMultivector, b) = index_to_bits(typeof(a), b)




#= MULTIVECTOR TYPE INFERENCE UTILITIES

The `<:AbstractMultivector` type system is complex -- any multivector type consists of four orthogonal aspects:

 - `multivectortype`: one of `Blade, Multivector, MixedMultivector`
 - `signature`: an all-bits type parameter defining the geometric algebra
 - `storagetype`: the container type `S` which stores components of a `<:CompositeMultivector{S}` (not applicable for `Blade`)
 - `eltype`: the numerical type of the components

The direction of promotion for these aspects are, respectively:

 - `Blade` -> `Multivector` -> `MixedMultivector`
 -  no promotion occurs -- an error is thrown if types of different signature are compared
 - `StaticVector` -> `Vector` -> `SparseVector` - the storage type shouldn't be affected by `Blade` types
 -  usual promotion behaviour of `<:Number` types

The master function `best_type` performs `<:AbstractMultivector` type promotion, and exposes ways of
controlling various type parameters (e.g., setting the grade or requiring a minimum `multivectortype` or `eltype`).
=#

multivectortype(::Type{<:Blade}) = Blade
multivectortype(::Type{<:Multivector}) = Multivector
multivectortype(::Type{<:MixedMultivector}) = MixedMultivector
multivectortype(::T) where {T<:AbstractMultivector} = multivectortype(T)

shared_sig(::Type{<:AbstractMultivector{sig}}...) where sig = sig
shared_sig(::Type{<:AbstractMultivector}...) = error("multivectors must share the same metric signature")
shared_sig(As::AbstractMultivector...) = shared_sig(typeof.(As)...)

set_eltype_parameter(::Type{<:Vector}, ::Type{T}) where {T} = Vector{T}
set_eltype_parameter(::Type{<:SparseVector{Tv,Ti} where Tv}, ::Type{T}) where {T,Ti} = SparseVector{T,Ti}
set_eltype_parameter(::Type{<:SVector}, ::Type{T}) where {T} = SVector{N,T} where N
set_eltype_parameter(::Type{<:MVector}, ::Type{T}) where {T} = MVector{N,T} where N

set_size_parameter(::Type{T}, ::Val{N}) where {T,N} = T
set_size_parameter(::Type{<:SVector{N′,T} where N′}, ::Val{N}) where {N,T} = SVector{N,T}
set_size_parameter(::Type{<:MVector{N′,T} where N′}, ::Val{N}) where {N,T} = MVector{N,T}


# `default_storagetype` should choose types appropriately by taking into account
# the algebra's dimension for optimal memory useage and performance
default_storagetype(::Type{Multivector}, T, dim) = dim >= 2^8 ? SparseVector{UInt,T} : Vector{T}
default_storagetype(::Type{MixedMultivector}, T, dim) = dim >= 8 ? SparseVector{UInt,T} : Vector{T}


# Ad hoc way of determining resulting storagetype from multiple types
# Algorithm: find the type which appears latest in this list...
const STORAGETYPES = [Nothing, StaticVector, Vector, SparseVector]
# ...and use that as the resulting storagetype.

storagetype(::Type{<:CompositeMultivector{S}}) where S = S
storagetype(::Type{<:Blade}) = Nothing
storagetype(::T) where {T<:AbstractMultivector} = storagetype(T)

function promote_storagetype(as)
	priority(a) = findfirst(S -> a <: S, STORAGETYPES)
	_, i = findmax(priority.(as))
	as[i]
end

function _best_type_parameters(::Type{MultivectorType}, a...;
		set_eltype::Type{SetT}=Nothing, promote_eltype_with::Type{PromT}=Union{}) where {MultivectorType,SetT,PromT}

	sig = shared_sig(a...)
	T = SetT == Nothing ? promote_type(eltype.(a)..., PromT) : SetT

	if MultivectorType == Blade
		return sig, T
	end

	S = promote_storagetype(storagetype.(a))
	if S === Nothing
		S = default_storagetype(MultivectorType, T, dimension(sig))
	end
	S = set_eltype_parameter(S, T)

	sig, T, S
end

_unwrap_type(::Type{Type{T}}) where T = T
_unwrap_type(T) = T

# memoized version for type stability
@generated function best_type_parameters(::Type{MultivectorType}, a...;
		set_eltype::Type{SetT}=Nothing, promote_eltype_with::Type{PromT}=Union{}) where {MultivectorType,SetT,PromT}
	_best_type_parameters(MultivectorType, _unwrap_type.(a)...; set_eltype=SetT, promote_eltype_with=PromT)
end

# best type is fastest when returning concrete types (as opposed to `UnionAll`s)
function best_type(::Type{Blade}, a...; grade::Val{k}=Val(missing), kwargs...) where k
	sig, T = best_type_parameters(Blade, a...; kwargs...)
	ismissing(k) ? Blade{sig,k,T} where {k} : Blade{sig,k,T}
end
function best_type(::Type{Multivector}, a...; grade::Val{k}=Val(missing), kwargs...) where k
	sig, T, S = best_type_parameters(Multivector, a...; kwargs...)
	if ismissing(k)
		Multivector{sig,k,S} where k
	else
		N = binomial(dimension(sig), k)
		Multivector{sig,k,set_size_parameter(S, Val(N))}
	end
end
function best_type(::Type{MixedMultivector}, a...; kwargs...)
	sig, T, S = best_type_parameters(MixedMultivector, a...; kwargs...)
	S = set_size_parameter(S, Val(2^dimension(sig)))
	MixedMultivector{sig,S}
end

# single argument form is used to modify parameters (e.g., to change the eltype)
best_type(a::Type{<:HomogeneousMultivector{sig,k}}; kwargs...) where {sig,k} = best_type(multivectortype(a), a; grade=Val(k), kwargs...)
best_type(a::Type{<:AbstractMultivector}; kwargs...) = best_type(multivectortype(a), a; kwargs...)
best_type(a::AbstractMultivector; kwargs...) = best_type(typeof(a); kwargs...)



#= CONVERSION

Conversion to more general types is possible:
 Blade ⊂ Multivector ⊂ MixedMultivector

The grade parameter of the target type is *not* respected. For subtypes of `HomogeneousMultivector{sig,k}`,
conversion to a type of a differing grade will still produce an object of grade `k`.
=#

# conversion of element/storage type
Base.convert(::Type{<:Blade{sig,k,T} where k}, a::Blade) where {sig,T} = Blade{sig,grade(a),T}(a.coeff, bitsof(a))
Base.convert(T::Type{<:Multivector}, a::Multivector) = best_type(T; grade=Val(grade(a)))(a.components)
Base.convert(::Type{<:MixedMultivector{sig,C}}, a::MixedMultivector) where {sig,C} = MixedMultivector{sig,C}(a.components)

# get (2^dim)-element vector of components for every basis blade in the algebra
full_components_vector(a::MixedMultivector{sig,<:AbstractVector}) where sig = a.components
function full_components_vector(a::Multivector)
	n, k = dimension(a), grade(a)
	fcv = zeros(eltype(a), 2^n) # TODO: use sparse vector? fcv is large and will remain mostly empty
	offset = multivector_index_offset(k, n)
	fcv[begin + offset : binomial(n, k) + offset] = a.components
	fcv
end
function full_components_vector(a::Blade)
	n = dimension(a)
	fcv = zeros(eltype(a), 2^n) # TODO: use sparse vector? fcv is large and will remain mostly empty
	fcv[bits_to_mmv_index(bitsof(a), n)] = a.coeff
	fcv
end



# conversion from lower multivector type
function Base.convert(T::Type{<:Multivector}, a::Blade)
	T = best_type(T; grade=Val(grade(a)))
	i = bits_to_index(T, bitsof(a))
	comps = [i == j ? convert(eltype(T), a.coeff) : zero(eltype(T)) for j ∈ 1:ncomponents(T)]
	T(comps)
end
Base.convert(T::Type{<:MixedMultivector}, a::AbstractMultivector) = T(full_components_vector(a))


# conversion from scalar
Base.convert(::Type{<:Blade{sig,k,T} where k}, a::Scalar) where {sig,T} = Blade{sig,0,T}(a, bits_scalar())
Base.convert(T::Type{<:Multivector}, a::Scalar) = best_type(T; grade=Val(0))([a])
Base.convert(T::Type{<:MixedMultivector}, a::Scalar) = zero(T) + a



#= PROMOTION

Promotion should result in objects of identical types *except* for the parameter {k}, which should not be changed.
Promotion should only occur between objects of the same signature -- otherwise, fall back on default promotion behaviour.

TODO: is it actually appropriate to define these promote rules? maybe it's more Julian to not do implicit conversions like this...
=#

# same multivector type
Base.promote_rule(T::Type{<:Blade{sig}}, S::Type{<:Blade{sig}}) where sig = best_type(Blade, T, S)
Base.promote_rule(T::Type{<:Multivector{sig}}, S::Type{<:Multivector{sig}}) where sig = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector{sig}}, S::Type{<:MixedMultivector{sig}}) where sig = best_type(MixedMultivector, T, S)

# promotion to higher multivector type
Base.promote_rule(T::Type{<:Multivector{sig}}, S::Type{<:Blade{sig}}) where sig = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector{sig}}, S::Type{<:AbstractMultivector{sig}}) where sig = best_type(MixedMultivector, T, S)

# promotion from scalars - is this a bad idea? (probably; not needed and explicit conversion is better)
# note that the grade parameter should *not* be preserved
# Base.promote_rule(T::Type{<:AbstractMultivector{sig}}, S::Type{<:Scalar}) where sig = best_type(multivectortype(T), T; promote_eltype_with=S)



#= CONVENIENCE CONSTRUCTORS AND CONVERTERS =#

blade_like(a::AbstractMultivector, coeff=1, bits::Unsigned=bits_scalar()) = Blade{signature(a),grade(bits),eltype(a)}(coeff, bits)

Multivector(a::Blade) = convert(best_type(Multivector, a; grade=Val(grade(a))), a)
MixedMultivector(a::HomogeneousMultivector) = convert(best_type(MixedMultivector, a), a)


"""
	mapcomponents(f, a::AbstractMultivector; kwargs...)

Apply a mapping `f :: Blade -> Scalar` component-wise to `a`.
The function `f` is given a `Blade` representing each component,
and should return a scalar (not a blade).

The same keyword arguments as [`best_type`](@ref) are accepted,
including `set_eltype` and `promote_eltype_with` which can be
used to specify the eltype of the resultant multivector.
"""
# mapcomponents(f, a::Blade; kwargs...) = best_type(a; kwargs...)(f(a)) # best_type preserves {k,bits}
mapcomponents(f, a::Blade; kwargs...) = best_type(a; kwargs...)(f(a), bitsof(a)) # best_type preserves {k,bits}
function mapcomponents(f, a::CompositeMultivector; kwargs...)
	a′ = zero(best_type(a; kwargs...))
	for b ∈ blades(a)
		a′.components[bits_to_index(a, bitsof(b))] = f(b)
	end
	a′
end


Base.float(T::Type{<:AbstractMultivector}) = best_type(T, set_eltype=float(eltype(T)))
Base.big(T::Type{<:AbstractMultivector}) = best_type(T, set_eltype=big(eltype(T)))
Base.complex(T::Type{<:AbstractMultivector}) = best_type(T, set_eltype=complex(eltype(T)))
Base.real(T::Type{<:AbstractMultivector}) = best_type(T, set_eltype=real(eltype(T)))

Base.float(a::AbstractMultivector) = convert(float(typeof(a)), a)
Base.big(a::AbstractMultivector) = convert(big(typeof(a)), a)
Base.complex(a::AbstractMultivector) = convert(complex(typeof(a)), a)

Base.real(a::AbstractMultivector) = mapcomponents(b -> real(b.coeff), a; set_eltype=real(eltype(a)))
Base.imag(a::AbstractMultivector) = mapcomponents(b -> imag(b.coeff), a; set_eltype=real(eltype(a)))



#= EQUALITY

Only equality between multivectors of the *almost the same type* is defined,
leveraging the promotion interface to enable comparison between differing types.

Multivector promotion (intentionally) does *not* unify the grade and bits parameters,
so the following methods must cover types which are identical except for those parameters.

Note that two blades of differing basis or two multivectors of differing grade are
treated as equal if they are both zero.
=#

==(a::Blade{sig}, b::Blade{sig}) where sig = bitsof(a) == bitsof(b) ? a.coeff == b.coeff : iszero(a) && iszero(b)
==(a::(Multivector{sig,k,S} where k), b::(Multivector{sig,k,S} where k)) where {sig,S} = grade(a) == grade(b) ? a.components == b.components : iszero(a) && iszero(b)
==(a::MixedMultivector{sig,S}, b::MixedMultivector{sig,S}) where {sig,S} = a.components == b.components

==(a::AbstractMultivector, b::Scalar) = iszero(b) && iszero(a) || isscalar(a) && scalar(a) == b
==(a::Scalar, b::AbstractMultivector) = iszero(a) && iszero(b) || isscalar(b) && a == scalar(b)

==(a::AbstractMultivector{sig}...) where sig = ==(promote(a...)...)
==(a::AbstractMultivector...) = false # multivectors with non-identical signatures are teated as non-equal


function Base.isapprox(a::Blade{sig}, b::Blade{sig}; kwargs...) where sig
	if bitsof(a) == bitsof(b)
		isapprox(a.coeff, b.coeff; kwargs...)
	else
		isapprox(a.coeff, zero(a.coeff); kwargs...) && isapprox(b.coeff, zero(a.coeff); kwargs...)
	end
end
Base.isapprox(a::T, b::T; kwargs...) where {T<:MixedMultivector} = isapprox(a.components, b.components; kwargs...)
function Base.isapprox(a::Multivector{sig}, b::Multivector{sig}; kwargs...) where sig
	if grade(a) == grade(b)
		isapprox(a.components, b.components; kwargs...)
	else
		# multivectors of different grade are approximately equal is they are both approximately zero
		isapprox(a, zero(a); kwargs...) && isapprox(b, zero(b); kwargs...)
	end
end
Base.isapprox(a::AbstractMultivector, b::AbstractMultivector; kwargs...) = isapprox(promote(a, b)...; kwargs...)

Base.isapprox(a::AbstractMultivector, b::Scalar; kwargs...) = isapprox(a, one(a)*b; kwargs...)
Base.isapprox(a::Scalar, b::AbstractMultivector; kwargs...) = isapprox(a*one(b), b; kwargs...)



#= COMPONENT ACCESS =#

blades(a::Blade) = (a,)
blades(a::Multivector{sig,k,<:AbstractVector}) where {sig,k} = (Blade{sig,k}(coeff, bits) for (coeff, bits) ∈ zip(a.components, bits_of_grade(k)))
blades(a::MixedMultivector{sig,<:AbstractVector}) where sig = (Blade{sig}(coeff, index_to_bits(a, i)) for (i, coeff) ∈ enumerate(a.components))
# blades(a::CompositeMultivector) = (Blade{signature(a)}(coeff, index_to_bits(a, i)) for (i, coeff) ∈ enumerate(a.components))
# blades(a::CompositeMultivector{<:AbstractDict}) where sig = (Blade{sig}(coeff, bits) for (bits, coeff) ∈ a.components)


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

function getcomponent(a::CompositeMultivector, I...)
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


