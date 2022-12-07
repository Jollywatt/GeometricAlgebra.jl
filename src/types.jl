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
- `T`: Numerical type of the scalar coefficient, retrieved with [`eltype()`](@ref).
"""
struct BasisBlade{Sig,K,T} <: AbstractMultivector{Sig}
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

# from scalar
BasisBlade{Sig}(coeff::T) where {Sig,T<:Scalar} = BasisBlade{Sig,0,T}(0, coeff)

# warning: doesn’t check that K == count_ones(bits)
BasisBlade{Sig,K}(pair::Pair) where {Sig,K} = let (bits, coeff) = pair
	BasisBlade{Sig,K,typeof(coeff)}(bits, coeff)
end



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
Multivector{Sig,K}(comps::S) where {Sig,K,S} = Multivector{Sig,K,S}(comps)



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
	ncomponents(::CompositeMultivector)

Number of independent components of a multivector instance (or type).
"""
ncomponents(a::Multivector) = length(a.comps)
ncomponents(a::Type{<:Multivector{Sig}}) where {Sig} = sum(ncomponents(dimension(a), k) for k in grade(a); init = 0)

Base.length(::AbstractMultivector) = error(
	"$length is not defined for multivectors. Do you mean $(repr(ncomponents))()?")


"""
	eltype(::AbstractMultivector)

The numerical type of the components of a multivector instance (or type).
"""
Base.eltype(::OrType{<:BasisBlade{Sig,K,T} where {Sig,K}}) where T = T
Base.eltype(::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = eltype(S)


"""
	grade(a)

Grade(s) of a multivector `a`.

The grade of a `BasisBlade{Sig,K}` or `Multivector{Sig,K}` is the second type parameter, `K`,
which may be an integer (if `a` is homogeneous) or a collection (a range or tuple of grades).

See also [`ishomogeneous`](@ref).
"""
grade(::OrType{<:BasisBlade{Sig,K}}) where {Sig,K} = K
grade(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = K

ishomogeneous(a) = isone(length(grade(a)))


bitsof(a::BasisBlade) = a.bits
componentbits(a::OrType{<:Multivector}) = componentbits(Val(dimension(a)), Val(grade(a)))

"""
	componentindex(a::Multivector, b)

Index of components vector of `a` corresponding to the `BasisBlade` or bits `b`.
"""
componentindex(a, b::BasisBlade) = componentindex(a, bitsof(b))
componentindex(a, bits) = findfirst(==(bits), componentbits(a))


function componentslice(a::OrType{<:Multivector}, k::Integer)
	i = findfirst(componentbits(a)) do bits
		count_ones(bits) == k
	end
	isnothing(i) && return 1:0
	range(i, length = ncomponents(dimension(a), k))
end


largest_type(::Multivector, ::BasisBlade) = Multivector
largest_type(::BasisBlade, ::Multivector) = Multivector
largest_type(::BasisBlade, ::BasisBlade) = BasisBlade





#= Grade inference =#

"""
	unify_grades(dim, p, q)

Return grade(s) containing both `p` and `q`.

In order to reduce the number of possible `Multivector` grade parameters each triggering
separate compilations, the result may be larger than the union `p ∪ q`.

# Examples
```jldoctest; setup = :(using GeometricAlgebra: unify_grades)
julia> unify_grades(4, 0:4, 2)
0:4

julia> unify_grades(4, 0, 2) # even multivectors are worth representing specifically
0:2:4

julia> unify_grades(4, (0, 4), 3) # not worth having a special type for grades (0, 3, 4)
0:4
```
"""
function unify_grades(dim, k)
	k = (0:dim) ∩ k
	length(k) == 1 ? first(k) : Tuple(k)
end
function unify_grades(dim, k::OrdinalRange)
	lo, hi = max(0, minimum(k)), min(dim, maximum(k))
	lo == hi && return lo
	Δ = abs(step(k))
	Δ == 1 ? (lo:hi) : (lo:Δ:hi)
end
function unify_grades(dim, p, q)
	p = unify_grades(dim, p)
	q = unify_grades(dim, q)

	p ⊆ q && return q
	p ⊇ q && return p
	
	p isa Integer && q isa Integer && minmax(p, q) == (0, dim) && return (0, dim)
	all(iseven, p) && all(iseven, q) && return 0:2:dim
	
	0:dim
end
unify_grades(dim, p, q, c...) = unify_grades(dim, unify_grades(dim, p, q), c...)


"""
	resulting_grades(combine, dim, p, q)

Non-zero grade(s) resulting from the application of `combine` on `dim`-dimensional multivectors of grade(s) `p` and `q`.
"""
resulting_grades(combine, dim, P, Q) = unify_grades(dim, (resulting_grades(combine, dim, p::Integer, q::Integer) for p in P, q in Q)...)





#= Constructors =#

constructor(::OrType{<:BasisBlade{Sig,K}}) where {Sig,K} = BasisBlade{Sig,K}
constructor(::OrType{<:Multivector{Sig,K}}) where {Sig,K} = Multivector{Sig,K}
Base.copy(a::Multivector) = constructor(a)(copy(a.comps))


function Base.similar(M::Type{Multivector{Sig,K}}, abc::OrType{<:AbstractMultivector}...) where {Sig,K}
	T = promote_type(eltype.(abc)...)
	C = componentstype(Sig, ncomponents(M), T)
	Multivector{Sig,K,C}
end

"""
	resulting_multivector_type(f, a, b, ...)

Multivector type with grade(s) and storage type appropriate for representing `f(a, b)`.
"""
function resulting_multivector_type(f, abc::OrType{<:AbstractMultivector{Sig}}...) where {Sig}
	dim = dimension(Sig)
	k = unify_grades(dim, resulting_grades(f, dim, grade.(abc)...))
	similar(Multivector{Sig,k}, abc...)
end


BasisBlade(a::BasisBlade) = a
Multivector(a::Multivector) = a
Multivector(a::BasisBlade{Sig,K}) where {Sig,K} = add!(zero(similar(Multivector{Sig,K}, a)), a)


Base.zero(::OrType{<:BasisBlade{Sig,K,T}}) where {Sig,K,T} = BasisBlade{Sig}(0 => numberzero(T))
Base.zero(a::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,K}(zeroslike(S, ncomponents(a)))
Base.one(::OrType{<:BasisBlade{Sig,K,T}}) where {Sig,K,T} = BasisBlade{Sig}(0 => numberone(T))
Base.one(::OrType{<:Multivector{Sig,K,S}}) where {Sig,K,S} = Multivector{Sig,0}(oneslike(S, 1))

Base.iszero(a::BasisBlade) = isnumberzero(a.coeff)
Base.iszero(a::Multivector) = all(isnumberzero, a.comps)
Base.isone(a::BasisBlade) = iszero(grade(a)) && isone(a.coeff)
Base.isone(a::Multivector) = isnumberone(a.comps[1]) && all(isnumberzero, a.comps[2:end])


Base.iseven(a::AbstractMultivector) = all(iseven, grade(a))
Base.isodd(a::AbstractMultivector) = all(isodd, grade(a))



#= Indexing and Iteration =#

"""
	a[k]
	grade(a::Multivector{Sig}, k) -> Multivector{Sig,k}

Grade `k` part of a multivetor `a`. Returns a zero `k`-vector if `k ∉ grade(a)`.

Multiple grades may be given with a range or tuple. The "grades" `+` and `-`
may be used as shortcuts for the even and odd parts, respectively.

# Examples
```jldoctest
julia> mv = Multivector{Cl(2,1), 0:4}(1:8)
8-component Multivector{Cl(2,1), 0:4, UnitRange{Int64}}:
 1
 2 v1 + 3 v2 + 4 v3
 5 v12 + 6 v13 + 7 v23
 8 v123

julia> grade(mv, 2)
3-component Multivector{Cl(2,1), 2, UnitRange{Int64}}:
 5 v12
 6 v13
 7 v23

julia> grade(mv, (0, 3))
2-component Multivector{Cl(2,1), (0, 3), Vector{Int64}}:
 1
 8 v123

julia> mv[+] # even part
4-component Multivector{Cl(2,1), 0:2:2, Vector{Int64}}:
 1
 5 v12 + 6 v13 + 7 v23
```

"""
grade(a::BasisBlade, k) = add!(zero(similar(Multivector{signature(a),k}, a)), a)

function grade(a::Multivector{Sig}, k::Integer) where {Sig}
	if k ∈ grade(a)
		Multivector{Sig,k}(view(a.comps, componentslice(a, k)))
	else
		zero(similar(Multivector{Sig,k}, a))
	end
end

function grade(a::Multivector{Sig}, K) where {Sig}
	K = unify_grades(dimension(Sig), K)
	Multivector{Sig,K}(vcat([grade(a, k).comps for k in K]...))
end

grade(a::Multivector, ::typeof(+)) = grade(a, 0:2:dimension(a))
grade(a::Multivector, ::typeof(-)) = grade(a, 1:2:dimension(a))

Base.getindex(a::AbstractMultivector, k) = grade(a, k)

eachgrade(a::Multivector) = (Multivector{signature(a),k}(view(a.comps, componentslice(a, k))) for k in grade(a))

scalar(a::BasisBlade{Sig,0}) where {Sig} = a.coeff
scalar(a::BasisBlade) = numberzero(eltype(a))
scalar(a::Multivector) = first(grade(a, 0).comps)

isscalar(a::Number) = true
isscalar(a::BasisBlade) = iszero(grade(a)) || isnumberzero(a)
function isscalar(a::Multivector)
	if 0 ∈ grade(a)
		if first(grade(a)) == 0
			all(isnumberzero, a.comps[2:end])
		else
			error("not implemented")
		end
	else
		iszero(a)
	end
end

blades(a::BasisBlade) = [a]
blades(a::Multivector{Sig}) where {Sig} = [BasisBlade{Sig}(bits => coeff) for (bits, coeff) ∈ zip(bitsof(a), a.comps)]

nonzero_components(a::BasisBlade) = isnumberzero(a.coeff) ? () : (bitsof(a) => a.coeff,)
nonzero_components(a::Multivector) = Iterators.filter(!isnumberzero∘last, zip(componentbits(a), a.comps))



function setindex!!(a::Multivector, val, i)
	if issetindexable(a.comps)
		setindex!(a.comps, val, i)
		a
	else
		constructor(a)(setindex(a.comps, val, i))
	end
end