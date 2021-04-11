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


"""
	HomogeneousMultivector{sig,k} <: AbstractMultivector{sig}

Supertype of grade `k` elements in the geometric algebra with metric signature `sig`.
"""
abstract type HomogeneousMultivector{sig,k} <: AbstractMultivector{sig} end

grade(::Type{<:HomogeneousMultivector{sig,k} where sig}) where k = k
grade(::T) where {T<:HomogeneousMultivector} = grade(T)


"""
	Blade{sig,k,bits,T} <: HomogeneousMultivector{sig,k}

A blade of grade `k` with basis blade `bits` and scalar coefficient of type `T`.

Parameters
===
- `sig`: metric signature defining the parent geometric algebra
- `k`: grade of the blade, equal to `grade(bits) === count_ones(bits)`
- `bits`: unsigned integer representing the indices of the basis vectors whose wedge product is the unit blade
- `T`: type of the scalar coefficient
"""
struct Blade{sig,k,bits,T} <: HomogeneousMultivector{sig,k}
	coeff::T
end
Blade{sig}(coeff::T, bits) where {sig,T} = Blade{sig,grade(bits),bits,T}(coeff)

# for internal use only (assumes `k == grade(bits)` without checking)
Blade{sig,k}(coeff::T, bits) where {sig,k,T} = Blade{sig,k,bits,T}(coeff)
Blade{sig,k,bits}(coeff::T) where {sig,k,bits,T} = Blade{sig,k,bits,T}(coeff)
# (Blade{sig,k,bits,T} where {k,bits})(coeff, bits) where {sig,T} = Blade{sig}(convert(T, coeff), bits)

"""
	bitsof(a::Blade{sig,k,bits,T}) = bits

The bits parameter of a blade, representing its basis `k`-blade.
"""
bitsof(::Type{<:Blade{sig,k,bits} where {sig,k}}) where bits = bits
bitsof(::T) where {T<:Blade} = bitsof(T)



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
const Scalar = Union{filter(T -> !(T <: AbstractMultivector), subtypes(Number))...}


"""
	eltype(a)

The numerical type of the components of a multivector or multivector type.
"""
Base.eltype(::Type{<:Blade{sig,k,bits,T} where {sig,k,bits}}) where T = T
Base.eltype(::Type{<:CompositeMultivector{S}}) where S = valtype(S)


ncomponents(::Type{<:Multivector{sig,k}}) where {sig,k} = binomial(dimension(sig), k)
ncomponents(::Type{<:MixedMultivector{sig}}) where sig = 2^dimension(sig)
ncomponents(::T) where {T<:AbstractMultivector} = ncomponents(T)



#= ZERO & ONE CONSTRUCTORS =#

zeroslike(::Type{Vector{T}}, I...) where T = zeros(T, I...)
zeroslike(::Type{SparseVector{Tv,Ti}}, I...) where {Tv,Ti} = spzeros(Tv, Ti, I...)
zeroslike(::Type{S}, I...) where {S<:StaticVector} = zeros(S)

Base.zero(::Type{<:Blade{sig,k,bits,T}}) where {sig,k,bits,T} = Blade{sig,k,bits,T}(zero(T))
Base.zero(::Type{<:Multivector{sig,k,S}}) where {sig,k,S} = Multivector{sig,k}(zeroslike(S, binomial(dimension(sig), k)))
Base.zero(::Type{<:MixedMultivector{sig,S}}) where {sig,S} = MixedMultivector{sig}(zeroslike(S, 2^dimension(sig)))

# needed to override default implementation for `<:Number`
Base.zero(::T) where {T<:AbstractMultivector} = zero(T)

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector{<:AbstractVector}) = iszero(a.components)
Base.iszero(a::CompositeMultivector{<:AbstractDict}) = error("unimplemented")

Base.one(a::Type{<:Blade}) = Blade{signature(a)}(one(eltype(a)), bits_scalar())



#= MULTIVECTOR TYPE INFERENCE UTILITIES

The `<:AbstractMultivector` type system is complex -- any multivector type consists of four orthogonal aspects:

 - `multivectortype`: one of `Blade, Multivector, MixedMultivector`
 - `signature`: an all-bits type parameter defining the geometric algebra
 - `storagetype`: the container type `S` which stores components of a `<:CompositeMultivector{S}`
 - `eltype`: the numerical type of the components

The direction of promotion for these aspects are, respectively:

 - `Blade` -> `Multivector` -> `MixedMultivector`
 -  no promotion occurs -- an error is thrown if types of different signature are compared
 - `Vector` -> `SparseVector`
 -  usual promotion of `<:Number` types

The master function `best_type` performs promotion, and exposes ways of
setting certain aspects of the resulting type (e.g., requiring a minimum
`multivectortype` or `eltype`).
=#

multivectortype(::Type{<:Blade}) = Blade
multivectortype(::Type{<:Multivector}) = Multivector
multivectortype(::Type{<:MixedMultivector}) = MixedMultivector
multivectortype(::T) where {T<:AbstractMultivector} = multivectortype(T)

shared_sig(::Type{<:AbstractMultivector{sig}}...) where sig = sig
shared_sig(::Type{<:AbstractMultivector}...) = error("multivectors must share the same metric signature")
shared_sig(T::AbstractMultivector...) = shared_sig(typeof.(T)...)


struct StorageType{S} end
unwrap(::Type{StorageType{S}}) where S = S

default_storagetype(sig, Tv, Ti) = Vector{Tv}
storagetype(::Type{<:Blade{sig,k,bits,T}}) where {sig,k,bits,T} = StorageType{default_storagetype(sig, T, typeof(bits))}
storagetype(::Type{<:Blade{sig,k,bits,T} where {k,bits}}) where {sig,T} = StorageType{default_storagetype(sig, T, UInt)}
storagetype(::Type{<:CompositeMultivector{S}}) where S = StorageType{S}
storagetype(::T) where {T<:AbstractMultivector} = storagetype(T)

Base.promote_rule(::Type{StorageType{Vector{T}}}, ::Type{StorageType{Vector{S}}}) where {T,S} = StorageType{Vector{promote_type(T, S)}}
Base.promote_rule(::Type{StorageType{SparseVector{Tv,Ti}}}, ::Type{StorageType{Vector{T}}}) where {Tv,Ti,T} = StorageType{SparseVector{promote_type(Tv, T),Ti}}
Base.promote_rule(::Type{StorageType{SparseVector{Tv,Ti}}}, ::Type{StorageType{SparseVector{Sv,Si}}}) where {Tv,Ti,Sv,Si} = StorageType{SparseVector{promote_type(Tv, Sv),promote_type(Ti, Si)}}
# fallback
Base.promote_rule(::Type{StorageType{T}}, ::Type{StorageType{S}}) where {T,S} = StorageType{promote_type(T, S)}

set_eltype_parameter(::Type{<:Vector}, ::Type{T}) where {T} = Vector{T}
set_eltype_parameter(::Type{<:SparseVector{Tv,Ti} where Tv}, ::Type{T}) where {T,Ti} = SparseVector{T,Ti}
set_eltype_parameter(::Type{<:SVector}, ::Type{T}) where {T} = SVector{N,T} where N
set_eltype_parameter(::Type{<:MVector}, ::Type{T}) where {T} = MVector{N,T} where N

set_size_parameter(::Type{T}, ::Val{N}) where {T,N} = T
set_size_parameter(::Type{<:SVector{N′,T} where N′}, ::Val{N}) where {N,T} = SVector{N,T}
set_size_parameter(::Type{<:MVector{N′,T} where N′}, ::Val{N}) where {N,T} = MVector{N,T}

# NOTE: best_type must work with non-concrete types in order for promotion to work with more than two objects.
# E.g., promote_type(x, x*y, x) == promote_type(promote_type(x, x*y), x)
#                                               \__ non-concrete __/
# Parameters which are allowed to be unspecified are {k,bits}. If any of {sig,S} are unspecified, do not use best_type to promote.

unwrap(::Type{Type{T}}) where T = T
unwrap(::Type{T}) where T = T

# must work where each `a...` is a multivector instance *or* a multivector type
@generated function best_type_parameters(a...; set_eltype::Type{SetT}=Nothing, promote_eltype_with::Type{PromT}=Union{}) where {SetT,PromT}
	a = unwrap.(a) # normalize Type{Type{A}} -> Type{A} 
	sig = shared_sig(a...)
	T = SetT == Nothing ? promote_type(eltype.(a)..., PromT) : SetT
	S = unwrap(promote_type(storagetype.(a)...))
	if S === StorageType error("storagetype promotion failed") end
	S = set_eltype_parameter(S, T)
	sig, T, S
end

# best type if fastest when returning concrete types (as opposed to `UnionAll`s)
function best_type(::Type{Blade}, a...; bits::Val{b}=Val(missing), kwargs...) where b
	sig, T = best_type_parameters(a...; kwargs...)
	ismissing(b) ? Blade{sig,k,bits,T} where {k,bits} : Blade{sig,grade(b),b,T}
end
function best_type(::Type{Multivector}, a...; grade::Val{k}=Val(missing), kwargs...) where k
	sig, T, S = best_type_parameters(a...; kwargs...)
	if ismissing(k)
		Multivector{sig,k,S} where k
	else
		N = binomial(dimension(sig), k)
		Multivector{sig,k,set_size_parameter(S, Val(N))}
	end
end
function best_type(::Type{MixedMultivector}, a...; kwargs...)
	sig, T, S = best_type_parameters(a...; kwargs...)
	S = set_size_parameter(S, Val(2^dimension(sig)))
	MixedMultivector{sig,S}
end

# single argument form is used to modify parameters (e.g., to change the eltype)
best_type(a::Type{<:Blade}; kwargs...) = best_type(Blade, a; bits=Val(bitsof(a)), kwargs...)
best_type(a::Type{<:Multivector}; kwargs...)  = best_type(Multivector, a; grade=Val(grade(a)), kwargs...)
best_type(a::Type{<:MixedMultivector}; kwargs...) = best_type(MixedMultivector, a; kwargs...)
best_type(a::AbstractMultivector; kwargs...) = best_type(typeof(a); kwargs...)




#= CONVERSION

Conversion to more general types is possible:
 Blade ⊂ Multivector ⊂ MixedMultivector

The grade parameter of the target type is *not* respected. For subtypes of `HomogeneousMultivector{sig,k}`,
conversion into a type of a differing grade will produce an object of grade `k`.
=#

# conversion of element/storage type
Base.convert(::Type{<:Blade{sig,k,bits,T} where {k,bits}}, a::Blade) where {sig,T} = Blade{sig,grade(a),bitsof(a),T}(a.coeff)
Base.convert(::Type{<:Multivector{sig,k,S} where k}, a::Multivector) where {sig,S} = add!(zero(Multivector{sig,grade(a),S}), a)
Base.convert(::Type{<:MixedMultivector{sig,C}}, a::MixedMultivector) where {sig,C} = add!(zero(MixedMultivector{sig,C}), a)

# conversion from lower multivector type
function Base.convert(T::Type{<:Multivector{sig,k,S} where k}, a::Blade) where {sig,S}
	T = Multivector{sig,grade(a),S}
	i = bits_to_key(T, bitsof(a))
	comps = [i == j ? convert(eltype(S), a.coeff) : zero(eltype(S)) for j ∈ 1:ncomponents(T)]
	T(comps)
end
function Base.convert(T::Type{<:MixedMultivector{sig,S}}, a::Blade) where {sig,S}
	i = 1 + bitsof(a)
	comps = [i == j ? convert(eltype(S), a.coeff) : zero(eltype(S)) for j ∈ 1:ncomponents(T)]
	MixedMultivector{sig,S}(comps)
end
Base.convert(T::Type{<:MixedMultivector}, a::HomogeneousMultivector) = add!(zero(T), a)

# conversion from scalar
Base.convert(::Type{<:Blade{sig,k,bits,T} where {k,bits}}, a::Scalar) where {sig,T} = Blade{sig,0,bits_scalar(),T}(a)
Base.convert(::Type{<:Multivector{sig,k,S} where k}, a::Scalar) where {sig,S} = zero(Multivector{sig,0,S}) + a
Base.convert(T::Type{<:MixedMultivector}, a::Scalar) = zero(T) + a



#= PROMOTION

Promotion should result in objects of identical types *except* for the parameters {k,bits}, which should not be changed.
Promotion should only occur between objects of the same signature -- otherwise, fall back on default promotion behaviour.
=#

# same multivector type
Base.promote_rule(T::Type{<:Blade{sig}}, S::Type{<:Blade{sig}}) where sig = best_type(Blade, T, S)
Base.promote_rule(T::Type{<:Multivector{sig}}, S::Type{<:Multivector{sig}}) where sig = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector{sig}}, S::Type{<:MixedMultivector{sig}}) where sig = best_type(MixedMultivector, T, S)

# promotion to higher multivector type
Base.promote_rule(T::Type{<:Multivector{sig}}, S::Type{<:Blade{sig}}) where sig = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector{sig}}, S::Type{<:AbstractMultivector{sig}}) where sig = best_type(MixedMultivector, T, S)

# promotion from scalars
# note that the grade/bits parameters should *not* be preserved
Base.promote_rule(T::Type{<:AbstractMultivector{sig}}, S::Type{<:Scalar}) where sig = best_type(multivectortype(T), T; promote_eltype_with=S)



#= CONVENIENCE CONSTRUCTORS AND CONVERTERS =#

blade_like(a::AbstractMultivector, coeff=1, bits::Unsigned=bits_scalar()) = Blade{signature(a),grade(bits),bits,eltype(a)}(coeff)

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
mapcomponents(f, a::Blade; kwargs...) = best_type(a; kwargs...)(f(a)) # best_type preserves {k,bits}
function mapcomponents(f, a::CompositeMultivector; kwargs...)
	a′ = zero(best_type(a; kwargs...))
	for b ∈ blades(a)
		a′.components[bits_to_key(a, bitsof(b))] = f(b)
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

Base.isapprox(a::AbstractMultivector, b::Scalar; kwargs...) = isapprox(promote(a, b)...; kwargs...)
Base.isapprox(a::Scalar, b::AbstractMultivector; kwargs...) = isapprox(promote(a, b)...; kwargs...)



#= COMPONENT ACCESS =#

blades(a::Blade) = (a,)
blades(a::Multivector{sig,k,<:AbstractVector}) where {sig,k} = (Blade{sig,k}(coeff, bits) for (coeff, bits) ∈ zip(a.components, FixedGradeBits(k)))
blades(a::MixedMultivector{sig,<:AbstractVector}) where sig = (Blade{sig}(coeff, unsigned(i - 1)) for (i, coeff) ∈ enumerate(a.components))
blades(a::CompositeMultivector{<:AbstractDict}) where sig = (Blade{sig}(coeff, bits) for (bits, coeff) ∈ a.components)


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

bits_to_key(::Type{<:Multivector{sig,k,<:AbstractVector}}, bits::Unsigned) where {sig,k} = bits_to_linear_index(bits)
bits_to_key(::Type{<:MixedMultivector{sig,<:AbstractVector}}, bits::Unsigned) where sig = 1 + bits
bits_to_key(::Type{<:CompositeMultivector{<:AbstractDict}}, bits::Unsigned) = bits
bits_to_key(::T, bits) where {T<:AbstractMultivector} = bits_to_key(T, bits)

parity_sign(I) = iseven(parity(sortperm(collect(I)))) ? +1 : -1

function getcomponent(a::Blade, I...)
	indices_to_bits(I) == bitsof(a) || return zero(eltype(a))
	s = parity_sign(I)
	s*a.coeff
end

getcomponent(a::Multivector, I...) = zero(eltype(a))
function getcomponent(a::Multivector{sig,k}, I::Vararg{Any,k}) where {sig,k}
	s = parity_sign(I)
	s*a.components[bits_to_key(a, indices_to_bits(I))]
end

function getcomponent(a::MixedMultivector, I...)
	s = parity_sign(I)
	s*a.components[bits_to_key(a, indices_to_bits(I))]
end

Base.getindex(a::AbstractMultivector, I::Integer...) = getcomponent(a, I...)
Base.getindex(a::AbstractMultivector) = getcomponent(a)

# experimental notations for grade selection
# Base.getindex(a::AbstractMultivector, I::Vector{<:Integer}) = grade(a, Iterators.only(I))
# Base.getindex(a::AbstractMultivector; grade) = GeometricAlgebra2.grade(a, grade)


