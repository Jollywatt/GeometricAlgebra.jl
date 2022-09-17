const OrType{T} = Union{T,Type{T}}

"""
	AbstractMultivector{Sig}

Supertype of all elements in the geometric algebra defined by the
metric signature `Sig`.

# Subtypes

```
                   AbstractMultivector{Sig}
                     /                  \\
   HomogeneousMultivector{Sig,K}    MixedMultivector{Sig,S}
       /                \\                             
Blade{Sig,K,T}    Multivector{Sig,K,S}                
                                                   
                  ╰───── CompositeMultivector{Sig,S} ─────╯
```

- `Blade`: a scalar multiple of a wedge product of orthogonal basis vectors.
- `Multivector`: a homogeneous multivector; a sum of same-grade blades.
- `MixedMultivector`: an inhomogeneous multivector. All elements in a geometric
   algebra can be represented as this type (though not most efficiently).

!!! note
	The mathematical definition of a ``k``-blade is the wedge product
	of ``k`` _vectors_, not necessarily basis vectors. Thus, not all
	``k``-blades are representable as a `Blade`, but are always representable
	as a sum of `Blade`s, or as a `Multivector`.

# Type Parameters

- `Sig`: The metric signature which defines the geometric algebra. This can be any
   all-bits value which satisfies the metric signature interface.
   For example, `(1, 1, 1)` or `EuclideanSignature(3)` both
   define the standard geometric algebra of ``ℝ^3``.
- `T`: The numerical type of the coefficient of a `Blade`.
- `K`: An `Int` specifying the grade of a `HomogeneousMultivector`.
- `S`: The storage type of the components of a `CompositeMultivector`. This is
   assumed to be mutable, and is usually a subtype of `Vector`, `MVector` or `SparseVector`.

"""
abstract type AbstractMultivector{Sig} end

Base.broadcastable(a::AbstractMultivector) = Ref(a)

"""
	signature(::AbstractMultivector{Sig}) -> Sig

The metric signature type parameter of the multivector object or type.
"""
signature(::OrType{<:AbstractMultivector{Sig}}) where {Sig} = Sig



"""
	HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig}

Abstract supertype of `Blade` and `Multivector`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade of the blade or multivector, retrieved with [`grade()`](@ref).
"""
abstract type HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig} end

"""
	grade(::HomogeneousMultivector{Sig,K}) -> K

The grade of a homogeneous multivector (`Blade` or `Multivector`) object or type.
"""
grade(::OrType{<:HomogeneousMultivector{Sig,K}}) where {Sig,K} = K



"""
	Blade{Sig,K,T} <: HomogeneousMultivector{Sig,K}

A blade of grade `K` with basis blade `bits` and scalar coefficient of type `T`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade of the blade, equal to `count_ones(bits)`, retrieved with [`grade()`](@ref).
- `T`: Numerical type of the scalar coefficient, retrieved with [`eltype()`](@ref).
"""
struct Blade{Sig,K,T} <: HomogeneousMultivector{Sig,K}
	bits::UInt
	coeff::T
end
Blade{Sig}(bits, coeff::T) where {Sig,T} = Blade{Sig,count_ones(bits),T}(bits, coeff)
Blade{Sig}(pair::Pair) where {Sig} = Blade{Sig}(pair...)

bitsof(a::Blade) = a.bits

mv_index(a::Blade) = bits_to_mv_index(bitsof(a))
mmv_index(a::Blade) = bits_to_mmv_index(bitsof(a), dimension(a))



"""
	Multivector{Sig,K,S} <: HomogeneousMultivector{Sig,K}

A homogeneous multivector of grade `K` with storage type `S`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade of the multivector, retrieved with [`grade()`](@ref).
- `S`: Storage type of the multivector components, usually a subtype of `AbstractVector`.
"""
struct Multivector{Sig,K,S} <: HomogeneousMultivector{Sig,K}
	components::S
end
Multivector{Sig,K}(comps::S) where {Sig,K,S} = Multivector{Sig,K,S}(comps)

mmv_offset(a::Multivector) = multivector_index_offset(grade(a), dimension(a))
mmv_slice(a::Multivector) = mmv_slice(grade(a), dimension(a))



"""
	MixedMultivector{Sig,S} <: AbstractMultivector{Sig}

A (possibly inhomogeneous) multivector.

All elements of a geometric algebra are representable as a `MixedMultivector`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `S`: Storage type of the multivector components, usually a subtype of `AbstractVector`.
"""
struct MixedMultivector{Sig,S} <: AbstractMultivector{Sig}
	components::S
end
MixedMultivector{Sig}(comps::S) where {Sig,S} = MixedMultivector{Sig,S}(comps)



const CompositeMultivector{Sig,S} = Union{Multivector{Sig,K,S},MixedMultivector{Sig,S}} where {K}



#= AbstractMultivector Interface =#


"""
	dimension(::AbstractMultivector)

The dimension of the underlying vector space of the geometric algebra.
See [`ncomponents`](@ref) for the dimension of the algebra (i.e., the
number of independent components of a general multivector).
"""
dimension(::AbstractMultivector{Sig}) where {Sig} = dimension(Sig)



"""
	ncomponents(::CompositeMultivector)

Number of independent components of a multivector object or type.

In ``n`` dimensions, this is ``\\binom{n}{k}`` for a `Multivector` and
``2^n`` for a `MixedMultivector`.
"""
ncomponents(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = ncomponents(Sig, K)
ncomponents(::OrType{<:MixedMultivector{Sig}}) where {Sig} = ncomponents(Sig)

Base.length(::AbstractMultivector) = error(
	"$length is not defined for multivectors. Do you mean $(repr(ncomponents))()?")

"""
	eltype(::AbstractMultivector)

The numerical type of the components of a multivector object or type.
"""
Base.eltype(::OrType{<:Blade{Sig,K,T} where {Sig,K}}) where T = T
Base.eltype(::OrType{<:CompositeMultivector{Sig,S}}) where {Sig,S} = eltype(S)


bitsof(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = bits_of_grade(K, dimension(Sig))
bitsof(::OrType{<:MixedMultivector{Sig}}) where {Sig} = mmv_index_to_bits(dimension(Sig))


largest_type(::MixedMultivector, ::AbstractMultivector) = MixedMultivector
largest_type(::Multivector, ::HomogeneousMultivector) = Multivector
largest_type(::Blade, ::Blade) = Blade
largest_type(a, b) = largest_type(b, a)






#= Constructors =#

constructor(::Blade{Sig,K}) where {Sig,K} = Blade{Sig,K}
constructor(::Multivector{Sig,K}) where {Sig,K} = Multivector{Sig,K}
constructor(::MixedMultivector{Sig}) where {Sig} = MixedMultivector{Sig}

Base.copy(a::CompositeMultivector) = constructor(a)(copy(a.components))

Base.zero(::OrType{<:Blade{Sig,K,T}}) where {Sig,K,T} = Blade{Sig}(0 => zero(T))
Base.zero(a::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,K}(zeroslike(S, ncomponents(a)))
Base.zero(a::OrType{<:MixedMultivector{Sig,S}}) where {Sig,S} = MixedMultivector{Sig}(zeroslike(S, ncomponents(a)))

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector) = all(iszero, a.components)

Base.one(::OrType{<:Blade{Sig,K,T}}) where {Sig,K,T} = Blade{Sig}(0 => one(T))
Base.one(::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,0}(oneslike(S, 1))
Base.one(a::OrType{<:MixedMultivector}) = zero(a) + one(eltype(a))

Base.isone(a::Blade) = iszero(grade(a)) && isone(a.coeff)
Base.isone(a::Multivector) = iszero(grade(a)) && isone(a.components[1])
Base.isone(a::MixedMultivector) = isone(a.components[1]) && all(iszero, a.components[2:end])



#= Conversion ==

Elements of a geometric algebra can be converted into 'larger' types
by calling the type constructor:

	Blade -> Multivector -> MixedMultivector

The eltype may be set as the second argument, as in `Multivector(a, Float32)`.
Converting to the same type returns identically, `Multivector(a::Multivector) === a`,
but a copy is *always* made when the second eltype argument is given.
=#

# Blade <- Blade
Blade(a::Blade) = a


# Multivector <- Blade
function Multivector(a::Blade{Sig,K,T′}, T=T′) where {Sig,K,T′}
	S = componentstype(Sig, ncomponents(Sig, K), T)
	if ismutabletype(S)
		add!(zero(Multivector{Sig,K,S}), a)
	else
		j = mv_index(a)
		comps = ntuple(i -> i == j ? T(a.coeff) : zero(T), length_from_type(S))
		Multivector{Sig,K,S}(S(comps))
	end
end

# Multivector <- Multivector
Multivector(a::Multivector) = a
function Multivector(a::Multivector{Sig,K}, T) where {Sig,K}
	S = with_eltype(typeof(a.components), T)
	if ismutabletype(S)
		add!(zero(Multivector{Sig,K,S}), a) # ensure a copy is made
	else
		Multivector{Sig,K,S}(convert(S, a.components))
	end
end


# MixedMultivector <- Blade
function MixedMultivector(a::Blade{Sig,K,T′}, T=T′) where {Sig,K,T′}
	S = componentstype(Sig, ncomponents(Sig), T)
	if ismutabletype(S)
		add!(zero(MixedMultivector{Sig,S}), a)
	else
		j = mmv_index(a)
		comps = ntuple(i -> i == j ? T(a.coeff) : zero(T), length_from_type(S))
		MixedMultivector{Sig,S}(S(comps))
	end
end

# MixedMultivector <- Multivector
function MixedMultivector(a::Multivector{Sig,K,S′}, T=eltype(S′)) where {Sig,K,S′}
	N = ncomponents(Sig)
	S = componentstype(Sig, N, T)
	if ismutabletype(S)
		add!(zero(MixedMultivector{Sig,S}), a) # ensure a copy is made
	else
		nbefore = mmv_offset(a) - 1
		nafter  = N - ncomponents(a) - nbefore

		comps = (ntuple(i -> zero(T), nbefore)..., a.components..., ntuple(i -> zero(T), nafter)...)
		MixedMultivector{Sig,S}(S(comps))
	end
end

# MixedMultivector <- MixedMultivector
MixedMultivector(a::MixedMultivector) = a
function MixedMultivector(a::MixedMultivector{Sig}, T) where {Sig}
	S = with_eltype(typeof(a.components), T)
	if ismutabletype(S)
		add!(zero(MixedMultivector{Sig,S}), a) # ensure a copy is made
	else
		MixedMultivector{Sig,S}(convert(S, a.components))
	end
end






#= Indexing and Iteration =#

grade(a::MixedMultivector, k) = Multivector{signature(a),k}(a.components[mmv_slice(k, dimension(a))])

scalarpart(a::Blade{Sig,0}) where {Sig} = a.coeff
scalarpart(a::Blade) = zero(eltype(a))
scalarpart(a::CompositeMultivector) = a.components[begin]

isscalar(a::Number) = true
isscalar(a::Blade{Sig,0}) where {Sig} = true
isscalar(a::Blade) = iszero(a)
isscalar(a::HomogeneousMultivector{Sig,0}) where {Sig} = true
isscalar(a::HomogeneousMultivector) = iszero(a)
isscalar(a::MixedMultivector) = all(iszero, a.components[2:end])


nonzero_components(a::Blade) = [bitsof(a) => a.coeff][iszero(a.coeff) ? [] : [1]]
nonzero_components(a::Multivector) = (bits => coeff for (bits, coeff) in zip(bits_of_grade(grade(a), dimension(a)), a.components) if !iszero(coeff))
nonzero_components(a::MixedMultivector) = (bits => coeff for (bits, coeff) in zip(mmv_index_to_bits(dimension(a)), a.components) if !iszero(coeff))