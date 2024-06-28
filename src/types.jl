const OrType{T} = Union{T,Type{T}}
const Scalar = Union{Number,SymbolicUtils.Symbolic}



"""
	AbstractMultivector{Sig}

Supertype of all elements in the geometric algebra defined by the
metric signature `Sig`.

# Subtypes

```
         AbstractMultivector{Sig}
            /               \\
BasisBlade{Sig,K,T}   Multivector{Sig,K,S}
```

- [`BasisBlade`](@ref): a scalar multiple of a wedge product of orthogonal basis vectors.
- [`Multivector`](@ref): a homogeneous or inhomogeneous multivector; a sum of basis blades.

"""
abstract type AbstractMultivector{Sig} end



"""
	BasisBlade{Sig,K,T}

A basis blade of grade `K` and scalar coefficient of type `T`.

Basis blades are scalar multiples of wedge products of orthogonal _basis_ vectors.

!!! note
	Not every ``k``-blade (i.e., wedge product of ``k`` linearly
	independent vectors) is representable as a `BasisBlade`.
	However, every ``k``-blade is a [`Multivector`](@ref) of grade `k`.

# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K::Int`: Grade of the blade, equal to `count_ones(bits)`, retrieved with [`grade()`](@ref).
- `T`: Numerical type of the scalar coefficient, retrieved with `eltype()`.
"""
struct BasisBlade{Sig,K,T} <: AbstractMultivector{Sig}
	coeff::T
	bits::UInt
end

"""
	BasisBlade{Sig}(bits, coeff)

Basis blade with indices encoded by `bits` and scalar coefficient `coeff`.

Indices are encoded in binary (e.g., ``v₁∧v₃∧v₄`` has bits `0b1101`).

# Examples
```jldoctest
julia> BasisBlade{3}(42, 0b110) # a grade 2 blade in 3 dimensions
BasisBlade{3, 2, Int64}:
 42 v23
```
"""
BasisBlade{Sig}(coeff::T, bits::Unsigned) where {Sig,T} = BasisBlade{Sig,count_ones(bits),T}(coeff, bits)

# warning: doesn’t check that K == count_ones(bits)
BasisBlade{Sig,K}(coeff::T, bits::Unsigned) where {Sig,K,T} = BasisBlade{Sig,K,T}(coeff, bits)

# from scalar
BasisBlade{Sig}(coeff::T) where {Sig,T<:Scalar} = BasisBlade{Sig,0,T}(coeff, UInt(0))

BasisBlade(a::BasisBlade) = a



"""
	Multivector{Sig,K,S} <: AbstractMultivector{Sig}

A general multivector with parts of grade `∈ K`.

For homogeneous `k`-vectors, the grade parameter `K` is an integer.
Inhomogeneous multivectors may be specified with a range or tuple of grades.


# Parameters
- `Sig`: Metric signature defining the geometric algebra, retrieved with [`signature()`](@ref).
- `K`: Grade(s) present in the multivector. Can be an integer or a collection of integers (a range or tuple).
- `S`: Storage type of the multivector components, usually a subtype of `AbstractVector`.
"""
struct Multivector{Sig,K,S} <: AbstractMultivector{Sig}
	comps::S
end

"""
	Multivector{Sig,K}(comps)

Multivector with grade(s) `K` and component vector `comps`.

Components are ordered first by grade, then lexicographically by bits (see [`componentbits`](@ref)).

# Examples
```jldoctest
julia> Multivector{3,0:3}(1:2^3)
8-component Multivector{3, 0:3, UnitRange{Int64}}:
 1
 2 v1 + 3 v2 + 4 v3
 5 v12 + 6 v13 + 7 v23
 8 v123

julia> grade(ans, 1)
3-component Multivector{3, 1, UnitRange{Int64}}:
 2 v1
 3 v2
 4 v3
```
"""
function Multivector{Sig,K}(comps::S) where {Sig,K,S}
	@assert length(comps) === ncomponents(Multivector{Sig,K}) """
	instances of $(Multivector{Sig,K}) have $(ncomponents(Multivector{Sig,K})) components, but received $(length(comps))"""
	Multivector{Sig,K,S}(comps)
end

function Multivector{Sig,K}(a::BasisBlade{Sig}) where {Sig,K}
	grade(a) ∈ K || error("$(constructor(a)) cannot be represented as a $(Multivector{Sig,K}), since $(grade(a)) ∉ $K")
	M = Multivector{Sig,K}
	i = componentindex(M, a)
	M(SingletonVector(a.coeff, i, ncomponents(M)))
end

Multivector(a::BasisBlade{Sig,K}) where {Sig,K} = Multivector{Sig,K}(a)
Multivector(a::Multivector) = a

Multivector{Sig,K}(a::Multivector{Sig,K}) where {Sig,K} = a
function Multivector{Sig,K′}(a::Multivector{Sig,K}) where {Sig,K,K′}
	K ⊆ K′ || error("$(constructor(a)) cannot be represented as a $(Multivector{Sig,K})")
	grade(a, K′)
end


Base.convert(::Type{Multivector{Sig,K,S}}, a::Multivector{Sig,K}) where {Sig,K,S} = Multivector{Sig,K}(convert(S, a.comps))
function Base.convert(::Type{Multivector{Sig,K,S}}, a::Multivector{Sig,K′}) where {Sig,K,K′,S}
	K′ ⊆ K || error("$(constructor(a)) cannot be represented as $(Multivector{Sig,K}), since $K′ ⊈ $K")
	add!(zero(Multivector{Sig,K,S}), a)
end
Base.convert(::Type{Multivector{Sig,K}}, a::Multivector{Sig,K′,S}) where {Sig,K,K′,S} = convert(Multivector{Sig,K,S}, a)



#= AbstractMultivector Interface =#

Base.broadcastable(a::AbstractMultivector) = Ref(a)


"""
	signature(::AbstractMultivector{Sig}) = Sig

The metric signature type parameter of the multivector instance (or type).
"""
signature(::OrType{<:AbstractMultivector{Sig}}) where {Sig} = Sig


"""
	dimension(sig)
	dimension(::AbstractMultivector)

The dimension of the underlying vector space of the geometric algebra.
See [`ncomponents`](@ref) for the dimension of the algebra (i.e., the
number of independent components of a general multivector).
"""
dimension(::OrType{<:AbstractMultivector{Sig}}) where {Sig} = dimension(Sig)


"""
	ncomponents(::Multivector)

Number of independent components of a multivector instance (or type).
"""
ncomponents(a::Multivector) = length(a.comps)
@static if VERSION >= v"1.8"
	Base.@assume_effects :foldable ncomponents(a::Type{<:Multivector}) = sum(ncomponents.(dimension(a), grade(a)); init = 0)
else
	ncomponents(a::Type{<:Multivector}) = sum(ncomponents.(dimension(a), grade(a)); init = 0)
end

Base.length(::AbstractMultivector) = error(
	"$length is not defined for multivectors. Do you mean $(repr(ncomponents))?")


Base.eltype(::OrType{<: BasisBlade{Sig,K,T} where {Sig,K}}) where {T} = T
Base.eltype(::OrType{<:Multivector{Sig,K,S} where {Sig,K}}) where {S} = eltype(S)


"""
	grade(a)

Grade or grades present in a multivector `a`.

The grade of a `BasisBlade{Sig,K}` or `Multivector{Sig,K}` is the second type parameter, `K`.
In the case of a multivector, `K` may be an integer (if it is homogeneous) or a collection (a range or tuple of grades).

See also [`ishomogeneous`](@ref).
"""
grade(::OrType{<: BasisBlade{Sig,K}}) where {Sig,K} = K
grade(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = K
grade(::Scalar) = 0

"""
	ishomogeneous(a)

Whether `a` is homogeneous, i.e., consists of nonzero parts of the same grade.
"""
ishomogeneous(a) = isone(length(grade(a)))

componentbits(a::OrType{<:Multivector}) = componentbits(Val(dimension(a)), Val(grade(a)))

"""
	componentindex(a::Multivector, b::Union{Unsigned,BasisBlade})

Index of the component of `a.comps` which corresponds to the basis blade `b`.
"""
componentindex(a::OrType{<:Multivector}, b::BasisBlade) = componentindex(a, b.bits)
componentindex(a::OrType{<:Multivector}, bits::Unsigned) = findfirst(==(bits), componentbits(a))


#= Constructors =#

constructor(::OrType{<:BasisBlade{Sig,K}}) where {Sig,K} = BasisBlade{Sig,K}
constructor(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = Multivector{Sig,K}
Base.copy(a::Multivector) = constructor(a)(copy(a.comps))


"""
	resulting_multivector_type(f, a, b, ...)

Return a `Multivector{Sig,K,S}` type with parameters (signature `Sig`, grade(s) `K` and storage type `S`)
appropriate for representing the result of `f(a, b)`.

Calls `resulting_grades(f, dimension(Sig), grade(a), grade(b), ...)` to determine `K`.
"""
function resulting_multivector_type(f, abc::OrType{<:AbstractMultivector{Sig}}...) where {Sig}
	dim = dimension(Sig)
	K = promote_grades(dim, resulting_grades(f, dim, grade.(abc)...))
	T = promote_type(eltype.(abc)...)
	S = componentstype(Sig, ncomponents(Multivector{Sig,K}), T)
	Multivector{Sig,K,S}
end





Base.zero(::OrType{<:BasisBlade{Sig,K,T} where K}) where {Sig,T} = BasisBlade{Sig}(numberzero(T))
Base.zero(a::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,K}(zeroslike(S, ncomponents(a)))

"""
	zero(::Type{Multivector{Sig,K,S})
	zero(::Type{Multivector{Sig,K}}, [T])

Multivector of metric signature `Sig` and grade(s) `K` with components all equal to zero.
If specified, the components array is of type `S`, or is the default array type with element type `T`.
"""
Base.zero(a::Type{Multivector{Sig,K}}, T::Type=Int) where {Sig,K} = zero(Multivector{Sig,K,componentstype(Sig, ncomponents(a), T)})

Base.one(::OrType{<:BasisBlade{Sig,K,T} where K}) where {Sig,T} = BasisBlade{Sig}(numberone(T))
Base.one(::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = add!(zero(Multivector{Sig,0 ∈ K ? K : 0,S}), numberone(eltype(S)), UInt(0))

Base.iszero(a::BasisBlade) = isnumberzero(a.coeff)
Base.iszero(a::Multivector) = all(isnumberzero, a.comps)
Base.isone(a::BasisBlade) = iszero(grade(a)) && isnumberone(a.coeff)
Base.isone(a::Multivector) = 0 ∈ grade(a) && isnumberone(a.comps[begin]) && all(isnumberzero, a.comps[begin + 1:end])

Base.iseven(a::AbstractMultivector) = all(iseven, grade(a))
Base.isodd(a::AbstractMultivector) = all(isodd, grade(a))

"""
	scalar(a) -> Number

The scalar component of a multivector.
"""
scalar(a::Scalar) = a
scalar(a::BasisBlade{Sig,0}) where {Sig} = a.coeff
scalar(a::BasisBlade) = numberzero(eltype(a))
scalar(a::Multivector) = 0 ∈ grade(a) ? a.comps[componentindex(a, UInt(0))] : numberzero(eltype(a))

"""
	isscalar(a)

Whether the only non-zero part of a multivector is its scalar part; `a == scalar(a)`.
"""
isscalar(a::Scalar) = true
isscalar(a::BasisBlade) = iszero(grade(a)) || isnumberzero(a)
function isscalar(a::Multivector)
	0 ∈ grade(a) || return all(isnumberzero, a.comps)
	for (coeff, bits) in nonzero_components(a)
		!iszero(bits) && return false
	end
	true
end




#= Indexing and Iteration =#


blades(a::BasisBlade) = Ref(a)
blades(a::Multivector{Sig}) where {Sig} = (BasisBlade{Sig}(coeff, bits) for (coeff, bits) ∈ zip(a.comps, componentbits(a)))

nonzero_components(a::BasisBlade) = isnumberzero(a.coeff) ? () : ((a.coeff, a.bits),)
nonzero_components(a::Multivector) = Iterators.filter(!isnumberzero∘first, zip(a.comps, componentbits(a)))
