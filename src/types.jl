const OrType{T} = Union{T,Type{T}}

"""
	AbstractMultivector{Sig}

Supertype of all elements in the geometric algebra defined by the
metric signature `Sig`.

# Subtypes

```
                        AbstractMultivector{Sig}
                          /                  \\
        HomogeneousMultivector{Sig,K}    Multivector{Sig,S}
            /               \\                             
BasisBlade{Sig,K,T}   KVector{Sig,K,S}                
                                                        
                      ╰─── CompositeMultivector{Sig,S} ───╯
```

- `BasisBlade`: a scalar multiple of a wedge product of orthogonal basis vectors.
- `KVector`: a homogeneous multivector; a sum of same-grade blades.
- `Multivector`: an inhomogeneous multivector. All elements in a geometric
   algebra can be represented as this type (though not always most efficiently).

# Type Parameters

- `Sig`: The metric signature which defines the geometric algebra. This can be any
   all-bits value which satisfies the metric signature interface.
   For example, `3` or `(1, 1, 1)` both define the standard geometric algebra over ``ℝ^3``.
- `T`: The numerical type of the coefficient of a `BasisBlade`.
- `K`: An `Int` specifying the grade of a `HomogeneousMultivector`.
- `S`: The storage type of the components of a `CompositeMultivector`, usually an `AbstractVector` subtype.

"""
abstract type AbstractMultivector{Sig} end

Base.broadcastable(a::AbstractMultivector) = Ref(a)

"""
	signature(::AbstractMultivector{Sig}) = Sig

The metric signature type parameter of the multivector instance (or type).
"""
signature(::OrType{<:AbstractMultivector{Sig}}) where {Sig} = Sig



"""
	HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig}

Abstract supertype of [`BasisBlade`](@ref) and [`KVector`](@ref).

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade of the blade or multivector, retrieved with [`grade()`](@ref).
"""
abstract type HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig} end

"""
	grade(::HomogeneousMultivector{Sig,K}) = K

The grade of a homogeneous multivector (a `BasisBlade` or `KVector`) instance (or type).
"""
grade(::OrType{<:HomogeneousMultivector{Sig,K}}) where {Sig,K} = K



"""
	BasisBlade{Sig,K,T} <: HomogeneousMultivector{Sig,K}

A basis blade of grade `K` and scalar coefficient of type `T`.

!!! note
	A `BasisBlade` represents a scalar multiple of a wedge product of orthogonal basis vectors.
	Not all ``k``-blades (in the sense of a wedge product of ``k`` linearly
	independent vectors) are representable as a `BasisBlade`. General ``k``-blades
	may be regarded as ``k``-vectors (see [`KVector`](@ref)).

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade of the blade, equal to `count_ones(bits)`, retrieved with [`grade()`](@ref).
- `T`: Numerical type of the scalar coefficient, retrieved with [`eltype()`](@ref).
"""
struct BasisBlade{Sig,K,T} <: HomogeneousMultivector{Sig,K}
	bits::UInt
	coeff::T
end
BasisBlade{Sig}(bits, coeff::T) where {Sig,T} = BasisBlade{Sig,count_ones(bits),T}(bits, coeff)
"""

	BasisBlade{Sig}(bits, coeff)
	BasisBlade{Sig}(bits => coeff)

Basis blade with indices encoded by `bits` and scalar coefficient `coeff`.

Indices are encoded in binary (e.g., ``v₁∧v₃∧v₄`` has bits `0b1101`).

# Examples
```jldoctest
julia> BasisBlade{3}(0b110 => 42) # a grade 2 blade in 3 dimensions
BasisBlade{3, 2, Int64}:
 42 v23
```
"""
BasisBlade{Sig}(pair::Pair) where {Sig} = BasisBlade{Sig}(pair...)

# warning: does’t check that K == count_ones(bits)
BasisBlade{Sig,K}(pair::Pair) where {Sig,K} = let (bits, coeff) = pair
	BasisBlade{Sig,K,typeof(coeff)}(bits, coeff)
end

bitsof(a::BasisBlade) = a.bits



"""
	KVector{Sig,K,S} <: HomogeneousMultivector{Sig,K}

A homogeneous multivector of grade `K` with storage type `S`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade of the multivector, retrieved with [`grade()`](@ref).
- `S`: Storage type of the multivector components, usually a subtype of `AbstractVector`.
"""
struct KVector{Sig,K,S} <: HomogeneousMultivector{Sig,K}
	comps::S
end

"""
	KVector{Sig,K}(comps)

KVector of grade `K` with components vector `comps`.

Components are ordered lexicographically by bits (see [`componentbits`](@ref)).

# Examples
```jldoctest
julia> KVector{3,2}(1:3) # 3D bivector
3-component KVector{3, 2, UnitRange{Int64}}:
 1 v12
 2 v13
 3 v23
```
"""
KVector{Sig,K}(comps::S) where {Sig,K,S} = KVector{Sig,K,S}(comps)
KVector(comps::AbstractVector) = KVector{length(comps),1}(comps)




"""
	Multivector{Sig,S} <: AbstractMultivector{Sig}

A general (possibly inhomogeneous) multivector.

All elements of a geometric algebra are representable as a `Multivector`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `S`: Storage type of the multivector components, usually a subtype of `AbstractVector`.
"""
struct Multivector{Sig,S} <: AbstractMultivector{Sig}
	comps::S
end

"""
	Multivector{Sig}(comps)

General multivector with components vector `comps`.

Components are ordered first by grade, then lexicographically by bits (see [`componentbits`](@ref)).

# Examples
```jldoctest
julia> Multivector{3}(1:2^3)
8-component Multivector{3, UnitRange{Int64}}:
 1
 2 v1 + 3 v2 + 4 v3
 5 v12 + 6 v13 + 7 v23
 8 v123

julia> grade(ans, 1)
3-component KVector{3, 1, UnitRange{Int64}}:
 2 v1
 3 v2
 4 v3
```
"""
Multivector{Sig}(comps::S) where {Sig,S} = Multivector{Sig,S}(comps)



const CompositeMultivector{Sig,S} = Union{KVector{Sig,K,S},Multivector{Sig,S}} where {K}

CompositeMultivector(a::BasisBlade) = KVector(a)
CompositeMultivector(a::CompositeMultivector) = a


#= AbstractMultivector Interface =#


"""
	dimension(sig)
	dimension(::AbstractMultivector)

The dimension of the underlying vector space of the geometric algebra.
See [`ncomponents`](@ref) for the dimension of the algebra (i.e., the
number of independent components of a general multivector).
"""
dimension(::OrType{<:AbstractMultivector{Sig}}) where {Sig} = dimension(Sig)



"""
	ncomponents(::CompositeMultivector)

Number of independent components of a multivector instance (or type).

In ``n`` dimensions, this is ``\\binom{n}{k}`` for a `KVector` and
``2^n`` for a `Multivector`.
"""
ncomponents(::OrType{<:KVector{Sig,K}}) where {Sig,K} = ncomponents(Sig, K)
ncomponents(::OrType{<:Multivector{Sig}}) where {Sig} = ncomponents(Sig)

Base.length(::AbstractMultivector) = error(
	"$length is not defined for multivectors. Do you mean $(repr(ncomponents))()?")

"""
	eltype(::AbstractMultivector)

The numerical type of the components of a multivector instance (or type).
"""
Base.eltype(::OrType{<:BasisBlade{Sig,K,T} where {Sig,K}}) where T = T
Base.eltype(::OrType{<:CompositeMultivector{Sig,S}}) where {Sig,S} = eltype(S)


bitsof(::OrType{<:KVector{Sig,K}}) where {Sig,K} = componentbits(Val(dimension(Sig)), Val(K))
bitsof(::OrType{<:Multivector{Sig}}) where {Sig} = componentbits(Val(dimension(Sig)))

bits_to_index(a::KVector, bits) = bits_to_kvector_index(bits)
bits_to_index(a::Multivector, bits) = bits_to_multivector_index(dimension(a), bits)

index_by_blade(a::KVector, b::BasisBlade) = bits_to_kvector_index(bitsof(b))
index_by_blade(a::Multivector{Sig}, b::BasisBlade{Sig}) where {Sig} = bits_to_multivector_index(bitsof(b), dimension(Sig))

"""
	componentindex(a::CompositeMultivector, b)

Index/indices of components vector of `a` corresponding to the `BasisBlade`/`KVector` `b`.
If `b` is a `BasisBlade` or `Unsigned` bits, returns the single index of the same component in `a`;
if `b` is a `KVector`, returns a `UnitRange` of the `grade(b)` components in `a`.
"""
componentindex(a, b::BasisBlade) = componentindex(a, bitsof(b))
componentindex(a::OrType{<:KVector}, bits::Unsigned) = bits_to_kvector_index(bits)
componentindex(a::OrType{<:Multivector}, bits::Unsigned) = bits_to_multivector_index(bits, dimension(a))
componentindex(a::OrType{<:Multivector{Sig}}, b::OrType{<:KVector{Sig}}) where {Sig} = multivector_slice(Val(dimension(Sig)), Val(grade(b)))


largest_type(::Multivector, ::AbstractMultivector) = Multivector
largest_type(::KVector, ::HomogeneousMultivector) = KVector
largest_type(::BasisBlade, ::BasisBlade) = BasisBlade
largest_type(a, b) = largest_type(b, a)






#= Constructors =#

constructor(::OrType{<:BasisBlade{Sig,K}}) where {Sig,K} = BasisBlade{Sig,K}
constructor(::OrType{<:KVector{Sig,K}}) where {Sig,K} = KVector{Sig,K}
constructor(::OrType{<:Multivector{Sig}}) where {Sig} = Multivector{Sig}

Base.copy(a::CompositeMultivector) = constructor(a)(copy(a.comps))

Base.zero(::OrType{<:BasisBlade{Sig,K,T}}) where {Sig,K,T} = BasisBlade{Sig}(0 => numberzero(T))
Base.zero(a::OrType{<:KVector{Sig,K,S}}) where {Sig,K,S} = KVector{Sig,K}(zeroslike(S, ncomponents(a)))
Base.zero(a::OrType{<:Multivector{Sig,S}}) where {Sig,S} = Multivector{Sig}(zeroslike(S, ncomponents(a)))

Base.iszero(a::BasisBlade) = isnumberzero(a.coeff)
Base.iszero(a::CompositeMultivector) = all(isnumberzero, a.comps)

Base.one(::OrType{<:BasisBlade{Sig,K,T}}) where {Sig,K,T} = BasisBlade{Sig}(0 => numberone(T))
Base.one(::OrType{<:KVector{Sig,K,S}}) where {Sig,K,S} = KVector{Sig,0}(oneslike(S, 1))
Base.one(a::OrType{<:Multivector}) = zero(a) + numberone(eltype(a))

Base.isone(a::BasisBlade) = iszero(grade(a)) && isone(a.coeff)
Base.isone(a::KVector) = iszero(grade(a)) && isone(a.comps[1])
Base.isone(a::Multivector) = isnumberone(a.comps[1]) && all(isnumberzero, a.comps[2:end])


#= Conversion ==

Elements of a geometric algebra can be converted into 'larger' types
by calling the type constructor:

	BasisBlade -> KVector -> Multivector

The eltype may be set as the second argument, as in `KVector(a, Float32)`.
Converting to the same type returns identically, `KVector(a::KVector) === a`,
but a copy is *always* made when the second eltype argument is given.
=#

# BasisBlade <- BasisBlade
BasisBlade(a::BasisBlade) = a


# KVector <- BasisBlade
KVector(a::BasisBlade) = KVector(a, numberorany(eltype(a)))
function KVector(a::BasisBlade{Sig,K}, T) where {Sig,K}
	n = ncomponents(Sig, K)
	S = componentstype(Sig, n, T)
	M = KVector{Sig,K,S}
	if issetindexable(S)
		add!(zero(M), a)
	else
		j = componentindex(M, a)
		M(S(ntuple(i -> i == j ? convert(T, a.coeff) : numberzero(T), n)))
	end
end

# KVector <- KVector
KVector(a::KVector) = a
function KVector(a::KVector{Sig,K}, T) where {Sig,K}
	S = componentstype(Sig, ncomponents(Sig, K), T)
	if issetindexable(S)
		add!(zero(KVector{Sig,K,S}), a) # ensure a copy is made
	else
		KVector{Sig,K,S}(convert(S, a.comps))
	end
end


# Multivector <- BasisBlade
Multivector(a::BasisBlade) = Multivector(a, numberorany(eltype(a)))
function Multivector(a::BasisBlade{Sig,K}, T) where {Sig,K}
	n = ncomponents(Sig)
	S = componentstype(Sig, n, T)
	M = Multivector{Sig,S}
	if issetindexable(S)
		add!(zero(M), a)
	else
		j = componentindex(M, a)
		M(S(ntuple(i -> i == j ? convert(T, a.coeff) : numberzero(T), n)))
	end
end

# Multivector <- KVector
Multivector(a::KVector) = Multivector(a, numberorany(eltype(a)))
function Multivector(a::KVector{Sig,K}, T) where {Sig,K}
	n = ncomponents(Sig)
	S = componentstype(Sig, n, T)
	if issetindexable(S)
		add!(zero(Multivector{Sig,S}), a) # ensure a copy is made
	else
		slice = multivector_slice(Val(dimension(Sig)), Val(K))
		nbefore = first(slice) - 1
		nafter = n - last(slice)

		comps = ntuple(i -> i ∈ slice ? a.comps[i - nbefore] : numberzero(T), Val(n))
		Multivector{Sig,S}(S(comps))
	end
end

# Multivector <- Multivector
Multivector(a::Multivector) = a
function Multivector(a::Multivector{Sig}, T) where {Sig}
	S = componentstype(Sig, ncomponents(Sig), T)
	if issetindexable(S)
		add!(zero(Multivector{Sig,S}), a) # ensure a copy is made
	else
		Multivector{Sig,S}(convert(S, a.comps))
	end
end



function Base.similar(M::Type{KVector{Sig,K}}, aa::AbstractMultivector...) where {Sig,K}
	T = promote_type(eltype.(aa)...)
	C = componentstype(Sig, ncomponents(M), T)
	KVector{Sig,K,C}
end

function Base.similar(M::Type{Multivector{Sig}}, aa::AbstractMultivector...) where {Sig}
	T = promote_type(eltype.(aa)...)
	C = componentstype(Sig, ncomponents(M), T)
	Multivector{Sig,C}
end


#= Indexing and Iteration =#

grade(a::BasisBlade{Sig,K}, k) where {Sig,K} = K == k ? a : zero(a)
grade(a::KVector{Sig,K,C}, k) where {Sig,K,C} = K == k ? a : zero(KVector{Sig,k,C})
grade(a::Multivector, k) = KVector{signature(a),k}(view(a.comps, multivector_slice(Val(dimension(a)), Val(k))))

scalar(a::BasisBlade{Sig,0}) where {Sig} = a.coeff
scalar(a::BasisBlade) = numberzero(eltype(a))
scalar(a::Multivector) = a.comps[begin]
scalar(a::KVector{Sig,0}) where {Sig} = a.comps[begin]
scalar(a::KVector{Sig}) where {Sig} = zero(eltype(a))

isscalar(a::Number) = true
isscalar(a::BasisBlade{Sig,0}) where {Sig} = true
isscalar(a::BasisBlade) = iszero(a)
isscalar(a::HomogeneousMultivector{Sig,0}) where {Sig} = true
isscalar(a::HomogeneousMultivector) = iszero(a)
isscalar(a::Multivector) = all(isnumberzero, a.comps[2:end])

blades(a::BasisBlade) = [a]
blades(a::CompositeMultivector{Sig}) where {Sig} = [BasisBlade{Sig}(bits => coeff) for (bits, coeff) ∈ zip(bitsof(a), a.comps)]

nonzero_components(a::BasisBlade) = isnumberzero(a.coeff) ? () : (bitsof(a) => a.coeff,)
nonzero_components(a::CompositeMultivector) = Iterators.filter(!isnumberzero∘last, zip(bitsof(a), a.comps))



function setindex!!(a::CompositeMultivector, val, i)
	if issetindexable(a.comps)
		setindex!(a.comps, val, i)
		a
	else
		constructor(a)(setindex(a.comps, val, i))
	end
end