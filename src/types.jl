const OrType{T} = Union{T,Type{T}}

"""
	AbstractMultivector{Sig}

Supertype of all elements in the geometric algebra defined by the
metric signature `Sig` (retrieved with the `signature` method).

Subtypes
--------

```
          AbstractMultivector
           /               \\
HomogeneousMultivector   MixedMultivector
   /       \\                         │
Blade   Multivector                  │
            │                        │
            ╰─ CompositeMultivector ─╯
```

- `Blade`: a scalar multiple of a wedge product of orthogonal basis vectors.
   Note that the mathematical definition of a ``k``-blade is the wedge product
   of ``k`` different _vectors_, not necessarily basis vectors. Thus, not all
   ``k``-blades are representable as a `Blade` (but always as a `Multivector`).
- `Multivector`: a homogeneous multivector; a sum of same-grade blades.
- `MixedMultivector`: an inhomogeneous multivector. All elements in a geometric
   algebra can be represented as this type.

"""
abstract type AbstractMultivector{Sig} end

Base.broadcastable(a::AbstractMultivector) = Ref(a)
Base.length(::AbstractMultivector) = error(
	"$length is not defined for multivectors. Do you mean $(repr(ncomponents))()?")


"""
	signature(::AbstractMultivector{Sig}) -> Sig

The metric signature type parameter of the multivector object or type.
"""
signature(::OrType{<:AbstractMultivector{Sig}}) where {Sig} = Sig


"""
	dimension(::AbstractMultivector)

The dimension of the underlying vector space of the geometric algebra.
See [`ncomponents`](@ref) the number of independent components of a multivector.
"""
dimension(::AbstractMultivector{Sig}) where {Sig} = dimension(Sig)

"""
	HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig}

Supertype of grade `K ∈ ℕ` elements in the geometric algebra with metric signature `Sig`.
"""
abstract type HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig} end

"""
	grade(::HomogeneousMultivector{Sig,K}) -> K

The grade of a homogeneous multivector (`Blade` or `Multivector`) object or type.
"""
grade(::OrType{<:HomogeneousMultivector{Sig,K}}) where {Sig,K} = K

"""
	Blade{Sig,K,T} <: HomogeneousMultivector{Sig,K}

A blade of grade `K ∈ ℕ` with basis blade `bits` and scalar coefficient of type `T`.

Parameters
----------
- `Sig`: metric signature defining the parent geometric algebra
- `K`: grade of the blade, equal to `count_ones(bits)`
- `T`: type of the scalar coefficient
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

A homogeneous multivector of grade `K ∈ ℕ` with storage type `S`.

Parameters
----------
- `Sig`: metric signature defining the parent geometric algebra
- `K`: grade of the multivector
- `S`: type in which the multivector components are stored; usually a vector-like or dictionary-like type
"""
struct Multivector{Sig,K,S} <: HomogeneousMultivector{Sig,K}
	components::S
end
Multivector{Sig,K}(comps::S) where {Sig,K,S} = Multivector{Sig,K,S}(comps)


"""
	MixedMultivector{Sig,S} <: AbstractMultivector{Sig}

A (possibly inhomogeneous) multivector.

All elements of a geometric algebra are representable as a `MixedMultivector`.

Parameters
----------
- `Sig`: metric signature defining the parent geometric algebra
- `S`: type in which the multivector components are stored; usually a vector-like or dictionary-like type
"""
struct MixedMultivector{Sig,S} <: AbstractMultivector{Sig}
	components::S
end
MixedMultivector{Sig}(comps::S) where {Sig,S} = MixedMultivector{Sig,S}(comps)


const CompositeMultivector{Sig,S} = Union{Multivector{Sig,K,S},MixedMultivector{Sig,S}} where {K}


#= AbstractMultivector Interface =#

mv_size(sig, k) = binomial(dimension(sig), k)
mmv_size(sig) = 2^dimension(sig)

"""
	ncomponents(::CompositeMultivector)

Number of independent components of a multivector object or type.

In ``n`` dimensions, this is ``\\binom{n}{k}`` for a `Multivector` and
``2^n`` for a `MixedMultivector`.
"""
ncomponents(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = mv_size(Sig, K)
ncomponents(::OrType{<:MixedMultivector{Sig}}) where {Sig} = mmv_size(Sig)

"""
	eltype(::AbstractMultivector)

The numerical type of the components of a multivector object or type.
"""
Base.eltype(::OrType{<:Blade{Sig,K,T} where {Sig,K}}) where T = T
Base.eltype(::OrType{<:CompositeMultivector{Sig,S}}) where {Sig,S} = eltype(S)


bitsof(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = bits_of_grade(K, dimension(Sig))
bitsof(::OrType{<:MixedMultivector{Sig}}) where {Sig} = mmv_index_to_bits(dimension(Sig))


#= Constructors =#

constructor(::Blade{Sig,K}) where {Sig,K} = Blade{Sig,K}
constructor(::Multivector{Sig,K}) where {Sig,K} = Multivector{Sig,K}
constructor(::MixedMultivector{Sig}) where {Sig} = MixedMultivector{Sig}

Base.zero(::OrType{<:Blade{Sig,K,T}}) where {Sig,K,T} = Blade{Sig}(0 => zero(T))
Base.zero(a::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,K}(zeroslike(S, ncomponents(a)))
Base.zero(a::OrType{<:MixedMultivector{Sig,S}}) where {Sig,S} = MixedMultivector{Sig}(zeroslike(S, ncomponents(a)))

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector) = iszero(a.components)

Base.one(::OrType{<:Blade{Sig,K,T}}) where {Sig,K,T} = Blade{Sig}(0 => one(T))
Base.one(::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,0}(oneslike(S, 1))
Base.one(a::OrType{<:MixedMultivector}) = add_scalar!(zero(a), 1)

Base.isone(a::Blade) = iszero(grade(a)) && isone(a.coeff)
Base.isone(a::Multivector) = iszero(grade(a)) && isone(a.components[])
Base.isone(a::MixedMultivector) = isone(a.components[1]) && iszero(a.components[2:end])

Base.copy(a::CompositeMultivector) = constructor(a)(copy(a.components))


#= Conversion ==

Multivectors can be converted into 'larger' types by calling the
type constructor:

	Blade -> Multivector -> MixedMultivector

The eltype may be set as the second argument, as in `Multivector(a, Float32)`.
Converting to the same type like `Multivector(a::Multivector)` returns `a`
identically, but a copy is always made when the eltype is given.
=#

largest_type(::MixedMultivector, ::AbstractMultivector) = MixedMultivector
largest_type(::Multivector, ::HomogeneousMultivector) = Multivector
largest_type(::Blade, ::Blade) = Blade
largest_type(a, b) = largest_type(b, a)



Blade(a::Blade) = a
Multivector(a::Multivector) = a
MixedMultivector(a::MixedMultivector) = a

function Multivector(a::Blade{Sig,K,T′}, T=T′) where {Sig,K,T′}
	C = componentstype(Sig, mv_size(Sig, K), T)
	add!(zero(Multivector{Sig,K,C}), a)
end
function Multivector(a::Multivector{Sig,K}, T) where {Sig,K}
	C = with_eltype(typeof(a.components), T)
	add!(zero(Multivector{Sig,K,C}), a)
end

function MixedMultivector(a::Blade{Sig,K,T′}, T=T′) where {Sig,K,T′}
	C = componentstype(Sig, mmv_size(Sig), T)
	add!(zero(MixedMultivector{Sig,C}), a)
end
function MixedMultivector(a::Multivector{Sig,K,C′}, T=eltype(C′)) where {Sig,K,C′}
	C = componentstype(Sig, mmv_size(Sig), T)
	add!(zero(MixedMultivector{Sig,C}), a)
end
function MixedMultivector(a::MixedMultivector{Sig}, T) where {Sig}
	C = with_eltype(typeof(a.components), T)
	add!(zero(MixedMultivector{Sig,C}), a)
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
isscalar(a::MixedMultivector) = iszero(a.components[2:end])

mmv_slice(a::Multivector) = mmv_slice(grade(a), dimension(a))

nonzero_components(a::Blade) = [bitsof(a) => a.coeff][iszero(a.coeff) ? [] : [1]]
nonzero_components(a::Multivector) = (bits => coeff for (bits, coeff) in zip(bits_of_grade(grade(a), dimension(a)), a.components) if !iszero(coeff))
nonzero_components(a::MixedMultivector) = (bits => coeff for (bits, coeff) in zip(mmv_index_to_bits(dimension(a)), a.components) if !iszero(coeff))