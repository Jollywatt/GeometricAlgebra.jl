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


# Ad hoc way of determining resulting storagetype from multiple types
# Algorithm: find the type which appears latest in this list...
const STORAGETYPES = [Nothing, StaticVector, Vector, SparseVector, Any]
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
	S === Nothing && (S = default_storagetype(sig, T))
	S = set_eltype_parameter(S, T)
	if ismissing(k)
		Multivector{sig,k,S} where k
	else
		N = binomial(dimension(sig), k)
		S = set_size_parameter(S, Val(N))
		Multivector{sig,k,S}
	end
end
function best_type(::Type{MixedMultivector}, a...; kwargs...)
	sig, T, S = best_type_parameters(MixedMultivector, a...; kwargs...)
	S === Nothing && (S = default_storagetype(sig, T))
	S = set_eltype_parameter(S, T)
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
function full_components_vector(a::Blade)
	n = dimension(a)
	fcv = zeros(eltype(a), 2^n) # TODO: use sparse vector? fcv is large and will remain mostly empty
	fcv[bits_to_mmv_index(bitsof(a), n)] = a.coeff
	fcv
end
function full_components_vector(a::Multivector)
	n, k = dimension(a), grade(a)
	fcv = zeros(eltype(a), 2^n) # TODO: use sparse vector? fcv is large and will remain mostly empty
	offset = multivector_index_offset(k, n)
	fcv[begin + offset : binomial(n, k) + offset] = a.components
	fcv
end
full_components_vector(a::MixedMultivector{sig,<:AbstractVector}) where sig = a.components



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
Base.convert(T::Type{<:MixedMultivector}, a::Scalar) = scalar_add!(zero(T), a)



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

