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
   ``k``-blades are representable as a `Blade`, but instead as a `Multivector`.
- `Multivector`: a homogeneous multivector; a sum of same-grade blades.
- `MixedMultivector`: an inhomogeneous multivector. All elements in a geometric
   algebra can be converted into this type.

=#

"""
	AbstractMultivector{sig}

Supertype of all elements in the geometric algebra over the vector space with
metric signature `sig`.
"""
abstract type AbstractMultivector{sig} <: Number end

dimension(::AbstractMultivector{sig}) where sig = dimension(sig)


"""
	HomogeneousMultivector{sig,k} <: AbstractMultivector{sig}

Supertype of grade `k` elements in the geometric algebra with metric signature `sig`.
"""
abstract type HomogeneousMultivector{sig,k} <: AbstractMultivector{sig} end

grade(::HomogeneousMultivector{sig,k}) where {sig,k} = k


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

# for internal use only (assumes `k == grade(bits)`)
Blade{sig,k}(coeff::T, bits) where {sig,k,T} = Blade{sig,k,bits,T}(coeff)
Blade{sig,k,bits}(coeff::T) where {sig,k,bits,T} = Blade{sig,k,bits,T}(coeff)

bitsof(::(Blade{sig,k,b} where {sig,k})) where b = b



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



const CompositeMultivector{S} = Union{Multivector{sig,k,S},MixedMultivector{sig,S}} where {sig,k}
const Scalar = Union{filter(T -> !(T <: AbstractMultivector), subtypes(Number))...}


"""
	eltype(a)

The numerical type of the components of a multivector or multivector type.
"""
Base.eltype(::Type{<:Blade{sig,k,bits,T} where {sig,k,bits}}) where T = T
Base.eltype(::Type{<:CompositeMultivector{S}}) where S = valtype(S)



#= ZERO & ONE CONSTRUCTORS =#

zeroslike(::Type{Vector{T}}, I...) where T = zeros(T, I...)
Base.zero(::Type{<:Blade{sig,k,bits,T}}) where {sig,k,bits,T} = Blade{sig,k,bits,T}(zero(T))
Base.zero(::Type{<:Multivector{sig,k,S}}) where {sig,k,S} = Multivector{sig,k}(zeroslike(S, binomial(dimension(sig), k)))
Base.zero(::Type{<:MixedMultivector{sig,S}}) where {sig,S} = MixedMultivector{sig}(zeroslike(S, 2^dimension(sig)))

# needed to override default implementation for `<:Number`
Base.zero(::T) where {T<:AbstractMultivector} = zero(T)

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector{<:AbstractVector}) = iszero(a.components)
Base.iszero(a::CompositeMultivector{<:AbstractDict}) = error("unimplemented")

Base.one(::Type{<:Blade{sig,k,bits,T}}) where {sig,k,bits,T} = Blade{sig,k,bits_scalar(),T}(one(T))


#= MULTIVECTOR TYPE INFERENCE UTILITIES =#

shared_sig(::AbstractMultivector{sig}...) where sig = sig
shared_sig(::AbstractMultivector...) = error("multivectors must share the same metric signature")
shared_sig(::Type{<:AbstractMultivector{sig}}...) where sig = sig

# function best_type(::Type{Blade}, a...; bits::Val{b}) where b
# 	sig = shared_sig(a...)
# 	T = promote_type(eltype.(a)...)
# 	Blade{sig,grade(b),b,T}
# end
# function best_type(::Type{Multivector}, a...; grade::Val{k}) where k
# 	sig = shared_sig(a...)
# 	T = promote_type(eltype.(a)...)
# 	Multivector{sig,k,Vector{T}}
# end
function best_type(::Type{Blade}, a...)
	sig = shared_sig(a...)
	T = promote_type(eltype.(a)...)
	Blade{sig,k,bits,T} where {k,bits}
end
function best_type(::Type{Multivector}, a...)
	sig = shared_sig(a...)
	T = promote_type(eltype.(a)...)
	Multivector{sig,k,Vector{T}} where k
end
function best_type(::Type{MixedMultivector}, a...)
	sig = shared_sig(a...)
	T = promote_type(eltype.(a)...)
	MixedMultivector{sig,Vector{T}}
end



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
Base.convert(T::Type{<:Multivector{sig,k,S} where k}, a::Blade) where {sig,S} = add!(zero(Multivector{sig,grade(a),S}), a)
Base.convert(T::Type{<:MixedMultivector}, a::HomogeneousMultivector) = add!(zero(T), a)

# conversion from scalar
Base.convert(::Type{<:Blade{sig,k,bits,T} where {k,bits}}, a::Scalar) where {sig,T} = Blade{sig,0,bits_scalar(),T}(a)
Base.convert(::Type{<:Multivector{sig,k,S} where k}, a::Scalar) where {sig,S} = zero(Multivector{sig,0,S}) + a
Base.convert(T::Type{<:MixedMultivector}, a::Scalar) = zero(T) + a


#= PROMOTION

Promotion should result in objects of identical types *except* for the grade parameter, which should not be changed.
=#

# same multivector type
Base.promote_rule(T::Type{<:Blade}, S::Type{<:Blade}) = best_type(Blade, T, S)
Base.promote_rule(T::Type{<:Multivector}, S::Type{<:Multivector}) = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector}, S::Type{<:MixedMultivector}) = best_type(MixedMultivector, T, S)

# promotion to higher multivector type
Base.promote_rule(T::Type{<:Multivector}, S::Type{<:Blade}) = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector}, S::Type{<:AbstractMultivector}) = best_type(MixedMultivector, T, S)

# promotion from scalars
Base.promote_rule(T::Type{<:AbstractMultivector}, S::Type{<:Scalar}) = error("unimplemented")

# hack: instances of `HomogeneousMultivector` which share identical types *expect* for a possibly
# differing grade or bits parameter should be treated as the "same" type as far as promotion is concerned
Base.promote(as::(Blade{sig,k,bits,T} where {k,bits})...) where {sig,T} = as
Base.promote(as::(Multivector{sig,k,S} where k)...) where {sig,S} = as








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

==(a::AbstractMultivector, b::Scalar) = isscalar(a) && scalar(a) == b
==(a::Scalar, b::AbstractMultivector) = isscalar(b) && a == scalar(b)


blades(a::Blade) = (a,)
blades(a::Multivector{sig,k,<:AbstractVector}) where {sig,k} = (Blade{sig,k}(coeff, bits) for (coeff, bits) ∈ zip(a.components, FixedGradeBits(k)))
blades(a::MixedMultivector{sig,<:AbstractVector}) where sig = (Blade{sig}(coeff, unsigned(i - 1)) for (i, coeff) ∈ enumerate(a.components))
blades(a::CompositeMultivector{<:AbstractDict}) where sig = (Blade{sig}(coeff, bits) for (bits, coeff) ∈ a.components)



