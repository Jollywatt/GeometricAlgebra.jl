#= MULTIVECTORS

There are three types representing objects in the geometric algebra, all
subtyping `AbstractMultivector`.
 - `Blade`: a scalar multiple of a wedge product of orthonormal basis vectors
 - `Multivector`: a homogeneous multivector; a sum of same-grade blades
 - `MixedMultivector`: an inhomogeneous multivector. All objects in the algebra
	can be expressed as a mixed multivector.
=#

"""
`AbstractMultivector{sig,C}`

Supertype of all elements in the geometric algebra in space of metric signature `sig`.

All multivector-type objects provide an association from unit blades ('ublades') to components.
The parameter `C`, which may be a vector or dict type, is the type in which the `ublade => comp`
association is represented.
The component type and ublade type can be determined with `eltype` and `keytype` respectively.

"""
abstract type AbstractMultivector{sig,C} end
# In earlier versions, AbstractMultivector <: Number. Until there's a case for why multivectors
# should subtype Number, I'll have them not be Numbers in order to avoid method
# ambiguities arising from many definitions like +(::AbstractMultivector, ::Number).

"""
```
Blade{sig,T,B,k} <: AbstractMultivector{sig,Pair{B,T}}
Blade{sig}(coeff, ublade)
```

A `k`-blade (i.e., wedge product of `k` basis vectors) with coefficient of type `T`
and unit blade represented by type `B`, in space of metric signature `sig`.

```
julia> Blade{(1,1,1)}(42, 0b101)
2-Blade{(1,1,1), Int64, UInt8, 2}
 42 v₁v₃

julia> Blade{(x=1,y=1,z=1)}(1, [:x, :z])
2-Blade{(x = 1, y = 1, z = 1), Int64, Vector{Symbol}, 2}
 1 xz
```
"""
struct Blade{sig,T,B,k} <: AbstractMultivector{sig,Pair{B,T}}
	coeff::T
	ublade::B
	function Blade{sig,T,B,k}(coeff, ublade) where {sig,T,B,k}
		@assert ublade_grade(ublade) == k
		new{sig,T,B,k}(coeff, ublade)
	end
end
Blade{sig}(coeff::T, ublade::B) where {sig,T,B} = Blade{sig,T,B,ublade_grade(ublade)}(coeff, ublade)


"""
```
Multivector{sig,C,k} <: AbstractMultivector{sig,C}
```

Homogeneous multivector of grade `k` in space of metric signature `sig`, containing
`binomial(dim(sig), k)` independent components.

Components are stored in an instance of type `C`, which may be an
`AbstractVector{T}` or an `AbstractDict{B,T}` where `T` is the component type.

When `C<:AbstractVector`, components are ordered lexicographically by lindex (linear index).
E.g., the component of a 3-multivector corresponding to `v1∧v3∧v4`
or `0b1101` is located at the `ublade2lindex(0b1101) == 3`rd index.
"""
struct Multivector{sig,C,k} <: AbstractMultivector{sig,C}
	comps::C
	function Multivector{sig,C,k}(comps) where {sig,C<:AbstractVector,k}
		@assert length(comps) == binomial(dim(sig), k)
		new{sig,C,k}(comps)
	end
	function Multivector{sig,C,k}(comps) where {sig,C<:AbstractDict,k}
		@assert all(ublade_grade(u) == k for u ∈ keys(comps))
		new{sig,C,k}(comps)
	end
end
Multivector{sig}(k::Integer, comps::C) where {sig,C} = Multivector{sig,C,k}(comps)

"""
```
Multivector{sig}(k, comps)
Multivector{sig}(comps::AbstractVector)
```

Construct a multivector of grade `k` from a vector or dict. When the argument `k`
is omitted, returns a 1-multivector (if `length(comps) == dim(sig)`).

Examples
===
```
julia> Multivector{(1,1,1)}(rand(3))
1-Multivector{(1, 1, 1), Vector{Float64}, 1}
 0.516268423028202 v₁
 0.280155164895841 v₂
 0.764458706918537 v₃

julia> Multivector{(x=1,y=1)}(2, [pi])
2-Multivector{(x = 1, y = 1), Vector{Irrational{:π}}, 2}
 π xy
```
"""
Multivector{sig}(comps::C) where {sig,C<:AbstractVector} = Multivector{sig,C,1}(comps)

"""
```
MixedMultivector{sig,C} <: AbstractMultivector{sig,C}
```

Inhomogeneous multivector in space of metric signature `sig`, containing `2^dim(sig)`
independent components.

Components are stored in an instance of type `C`, which may be an
`AbstractVector{T}` or an `AbstractDict{B,T}` where `T` is the component type,
similar to `Multivectors{sig,C,k}`.

When `C<:AbstractVector`, components are ordered by the binary value of the unit blade.
E.g., the component corresponding to `v1∧v3∧v4` is at the `Int(0b1101) == 13`th index.
"""
struct MixedMultivector{sig,C} <: AbstractMultivector{sig,C}
	comps::C
	function MixedMultivector{sig,C}(comps) where {sig,C<:AbstractVector}
		@assert length(comps) == 2^dim(sig)
		new{sig,C}(comps)
	end
	function MixedMultivector{sig,C}(comps) where {sig,C<:AbstractDict}
		new{sig,C}(comps)
	end
end
MixedMultivector{sig}(comps::C) where {sig,C} = MixedMultivector{sig,C}(comps)


const HomogeneousMultivector{k} = Union{Blade{sig,T,B,k},Multivector{sig,C,k}} where {sig,T,B,C}
const CompositeMultivector{sig,C} = Union{Multivector{sig,C},MixedMultivector{sig,C}}
const multivector_types = [Blade, Multivector, MixedMultivector]



#TODO: equality. Should == involve implicit conversions?
for eq ∈ [:(==), :(≈)]
	@eval $eq(a::Blade{sig}, b::Blade{sig}; kwargs...) where sig = a.ublade == b.ublade && $eq(a.coeff, b.coeff; kwargs...)
	@eval $eq(a::T, b::T; kwargs...) where T<:CompositeMultivector{sig,C} where {sig,C<:AbstractVector} = $eq(a.comps, b.comps; kwargs...)
	@eval function $eq(a::T, b::T; kwargs...) where T<:CompositeMultivector{sig,C} where {sig,C<:AbstractDict}
		for u ∈ union(keys(a.comps), keys(b.comps))
			$eq(getcomp(a, u), getcomp(b, u); kwargs...) ? continue : return false
		end
		true
	end
end

"""
`signature(type)`

The metric signature of multivectors of type `type`. The definition
`signature(x) = signature(typeof(x))` is provided so that instances may be passed.
"""
signature(::Type{<:AbstractMultivector{sig}}) where sig = sig
signature(a::AbstractMultivector) = signature(typeof(a))

# MULTIVECTOR TYPE GYMNASTICS

# assert that multivectors share same metric signature and return it
shared_sig(::Type{<:AbstractMultivector{sig}}...) where sig = sig
shared_sig(::Type{<:AbstractMultivector}...) = error("multivectors must share the same metric signature")
shared_sig(as::AbstractMultivector...) = shared_sig(typeof.(as)...)

# give the most sensible storage type of a CompositeMultivector with given signature, eltype, keytype
storagetype(sig, T, B::Type{<:Unsigned}) = dim(sig) <= 8 ? Vector{T} : SparseVector{T}
storagetype(sig, T, B) = Dict{B,T}

# this might be simplified if BitVectors are banished
promote_ublade_type(as::Type{<:Unsigned}...) = promote_type(as...)
for T ∈ [Vector{<:Integer}, BitVector]
	@eval promote_ublade_type(a::Type{<:Unsigned}, b::Type{<:$T}) = $T
	@eval promote_ublade_type(a::Type{<:$T}, b::Type{<:Unsigned}) = $T
	@eval promote_ublade_type(a::Type{<:$T}, b::Type{<:$T}) = $T
end

promote_ublade_type(::Type{T}...) where T = T

# give the most sensible CompositeMultivector type which can represent the given arguments
function best_type(M::Type{Multivector}, as::Type{<:AbstractMultivector}...; k)
	sig = shared_sig(as...)
	# T = promote_type(filter(t -> !(t<:Derivative), eltype.(as))...)
	T = promote_type(eltype.(as)...)
	B = promote_ublade_type(keytype.(as)...)
	C = storagetype(sig, T, B)
	Multivector{sig,C,k}
end
function best_type(M::Type{MixedMultivector}, as::Type{<:AbstractMultivector}...)
	sig = shared_sig(as...)
	# T = promote_type(filter(t -> !(t<:Derivative), eltype.(as))...)
	T = promote_type(eltype.(as)...)
	B = promote_ublade_type(keytype.(as)...)
	C = storagetype(sig, T, B)
	MixedMultivector{sig,C}
end
# function best_type(M, a::AbstractMultivector{})
best_type(M, as::AbstractMultivector...; args...) = best_type(M, typeof.(as)...; args...)


function Base.convert(::Type{MixedMultivector}, a::AbstractMultivector)
	b = zero(best_type(MixedMultivector, a))
	for (u, v) ∈ comps(a)
		setcomp!(b, u, v)
	end
	b
end

function Base.convert(::Type{Multivector}, a::Union{Blade,Multivector})
	k = grade(a)
	b = zero(best_type(Multivector, a; k))
	for (u, v) ∈ comps(a)
		setcomp!(b, u, v)
	end
	b
end




best_ublade_type(n) = n <= 8 ? UInt8 : n <= 16 ? UInt16 : n <= 32 ? UInt32 : n <= 64 ? UInt64 : n <= 128 ? UInt128 : Vector{Int}

Base.eltype(::Type{<:Blade{sig,T}}) where {sig,T} = T
Base.eltype(::Type{<:CompositeMultivector{sig,<:AbstractVector{T}}}) where {sig,T} = T
Base.eltype(::Type{<:CompositeMultivector{sig,<:AbstractDict{B,T}}}) where {sig,B,T} = T

Base.keytype(::Type{<:Blade{sig,T,B}}) where {sig,T,B} = B
Base.keytype(::Type{<:CompositeMultivector{sig,<:AbstractVector}}) where {sig} = best_ublade_type(dim(sig))
Base.keytype(::Type{<:CompositeMultivector{sig,<:AbstractDict{B}}}) where {sig,B} = B
Base.keytype(a::AbstractMultivector) = keytype(typeof(a))


"""
```
basis(::Type{<:AbstractMultivector}, i)
basis(::Type{<:AbstractMultivector})
```

Return the `i`th basis element of the given type, or a tuple of all basis elements if no
index is given.
For homogeneous types with grade `k`, this gives the `i`th basis `k`-blade in lexicographic order.
For mixed multivector types, this gives the basis vectors / 1-blades.

Examples
===
```
julia> x, y, z = basis(Blade{(1,1,1),Float64,UInt,1})
3-element Vector{Blade{(1, 1, 1), Float64, UInt64, 1}}:
 1.0 v₁
 1.0 v₂
 1.0 v₃

julia> basis(Blade{(1,1,1),Float64,UInt,1}, 2) == y
true

julia> basis(Multivector{(1,1,1),Vector{Float},2}, 1) == x*y
true
```

"""
basis(::Type{Blade{sig,T,B,1}}, i) where {sig,T,B} = Blade{sig,T,B,1}(one(T), convert_ublade(B, [i]))
basis(::Type{Blade{sig,T,B,k}}, i) where {sig,T,B,k} = Blade{sig,T,B,k}(one(T), convert_ublade(B, lindex2ublade(k, i)))
function basis(M::Type{Multivector{sig,C,k}}, i) where {sig,C<:AbstractVector,k}
	a = zero(M)
	setcomp!(a, lindex2ublade(k, i), one(eltype(M)))
	a
end
basis(::Type{Multivector{sig,C}}, i) where {sig,C<:AbstractVector} = basis(Multivector{sig,C,1}, i)
function basis(M::Type{MixedMultivector{sig,C}}, i) where {sig,C<:AbstractVector}
	a = zero(M)
	setcomp!(a, ublade_bv(sig, keytype(M), i), one(eltype(M)))
	a
end
basis(M::Type{<:CompositeMultivector{sig,C}}, i) where {sig,C<:AbstractDict} = M(C(ublade_bv(sig, keytype(M), i) => one(eltype(M))))

basis(M::Type{<:AbstractMultivector{sig}}) where sig = [basis(M, i) for i ∈ 1:dim(sig)]
basis(M::Type{<:HomogeneousMultivector{k}}) where k = [basis(M, i) for i ∈ 1:binomial(dim(signature(M)), k)]

# basis(sig::Type{<:AbstractMetricSignature}, i) = basis(sig(), i)
"""
```
basis(sig, i)
basis(sig)
```

Return a basis vector `i` for the space of metric signature `sig`,
or a tuple of all basis vectors (if `sig` has specified dimension).

Examples
===
```
julia> basis((x=1, y=1))
2-element Vector{Blade{(x = 1, y = 1), Float64, UInt64, 1}}:
 1.0 x
 1.0 y

julia> basis((-1,1,1,1), 2)
1-Blade{(-1, 1, 1, 1), Float64, UInt64, 1}
 1.0 v₂

julia> basis(EuclideanSignature, :t)
1-Blade{EuclideanSignature, Float64, Vector{Symbol}, 1}
 1.0 t
```
"""
basis(sig, i) = basis(Blade{sig,Float64,UInt,1}, i)
basis(sig) = [basis(sig, i) for i ∈ 1:dim(sig)]

# good? bad? generate all 2^dim(sig) basis multivectors?
fullbasis(sig) = (Blade{sig}(1, UInt(i - 1)) for i ∈ 1:2^dim(sig))
basis(sig, i::Symbol) = Blade{sig,Float64,Vector{Symbol},1}(1, [i])

"""
```
vol(x)
vol(::Type)
`

The volume form or psuedoscalar element.
"""
# @pure ??
vol(T::Type{<:AbstractMultivector{sig}}) where sig = Blade{sig}(one(eltype(T)), ublade_vol(sig, keytype(T)))
vol(::T) where T<:AbstractMultivector = vol(T)

# COMPONENT ACCESS

getcomp(a::Blade{sig,T}, ublade) where {sig,T} = a.ublade == ublade ? a.coeff : zero(T)

function getcomp(a::Multivector{sig,<:AbstractVector,k}, ublade::Unsigned) where {sig,k}
	ublade_grade(ublade) == k || return zero(eltype(a))
	i = ublade2lindex(ublade)
	a.comps[i]
end

getcomp(a::MixedMultivector{sig,<:AbstractVector}, ublade::Unsigned) where sig = a.comps[begin + ublade]
# this is a headache with the zoo of ublade types
getcomp(a::CompositeMultivector{sig,<:AbstractDict}, ublade) where sig = get(a.comps, ublade, zero(eltype(a)))


function setcomp!(a::Multivector{sig,<:AbstractVector,k}, ublade, v) where {sig,k}
	ublade_grade(ublade) == k || error("cannot set non-$k grade components of $k-multivector")
	i = ublade2lindex(ublade)
	a.comps[i] = v
end

setcomp!(a::MixedMultivector{sig,<:AbstractVector}, ublade, v) where sig = a.comps[begin + ublade] = v
setcomp!(a::CompositeMultivector{sig,<:AbstractDict}, ublade, v) where sig = a.comps[ublade] = v
addcomp!(a::AbstractMultivector, ublade, v) = setcomp!(a, ublade, getcomp(a, ublade) + v)


"""
`comps(a::AbstractMultivector)`

Iterate components of multivector, generating `ublade => coeff` pairs (even when
coefficients are stored in a vector).
"""
comps(a::Blade) = (a.ublade => a.coeff,)

#TODO: this could be made more efficient with an iterator which calls next_bit_permutation once each time
comps(a::Multivector{sig,C,k}) where {sig,C<:AbstractVector,k} = (lindex2ublade(k, i) => coeff for (i, coeff) ∈ enumerate(a.comps))

comps(a::MixedMultivector{sig,C}) where {sig,C<:AbstractVector} = (unsigned(i - 1) => coeff for (i, coeff) ∈ enumerate(a.comps))
comps(a::CompositeMultivector{sig,C}) where {sig,C<:AbstractDict} = a.comps



# map components of multivector without mixing ublades
function mapcomp!(f::Function, a::CompositeMultivector)
	for (u, v) ∈ comps(a)
		setcomp!(a, u, f(u, v))
	end
	a
end
function mapcomp(f::Function, a::Blade)
	Blade{signature(a)}(f(a.ublade, a.coeff), a.ublade)
end
mapcomp(f::Function, a::AbstractMultivector) = mapcomp!(f, deepcopy(a))



_zeros(::Type{Vector{T}}, n) where T = zeros(T, n)
_zeros(::Type{SparseVector{T}}, n) where T = spzeros(T, n)

Base.zero(::Type{<:Blade{sig,T,B,k}}) where {sig,T,B,k} = Blade{sig}(zero(T), ublade_scalar(sig, B))
Base.zero(::Type{<:Multivector{sig,C,k}}) where {sig,C<:AbstractVector,k} = Multivector{sig}(k, _zeros(C, binomial(dim(sig), k)))
Base.zero(::Type{<:MixedMultivector{sig,C}}) where {sig,C<:AbstractVector} = MixedMultivector{sig}(_zeros(C, 2^dim(sig)))
Base.zero(T::Type{<:CompositeMultivector{sig,C}}) where {sig,C<:AbstractDict} = T(C())
Base.zero(a::AbstractMultivector) = zero(typeof(a)) # not needed if AbstractMultivector <: Number

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector{sig,C}) where {sig,C<:AbstractVector} = iszero(a.comps)
Base.iszero(a::CompositeMultivector{sig,C}) where {sig,C<:AbstractDict} = all(iszero(v) for (u, v) ∈ a.comps)


Base.one(::Type{<:Blade{sig,T,B}}) where {sig,T,B} = Blade{sig}(one(T), ublade_scalar(sig, B))
Base.oneunit(a::Blade{sig,T,B,k}) where {sig,T,B,k} = Blade{sig}(one(T), a.ublade)
Base.one(a::AbstractMultivector) = one(typeof(a)) # not needed if AbstractMultivector <: Number

function Base.one(T::Type{<:Multivector{sig,C}}) where {sig,C}
	a = zero(Multivector{sig,C,0}) # scalar / 0-grade
	setcomp!(a, ublade_scalar(sig, keytype(T)), one(eltype(T)))
	a
end

function Base.one(T::Type{<:MixedMultivector})
	a = zero(T)
	setcomp!(a, ublade_scalar(signature(T), keytype(T)), one(eltype(T)))
	a
end


# GRADE SELECTION

"""
`grade(a::AbstractMultivector)`

The grade of the multivector `a`. For `MixedMultivector`s, returns a vector of each non-zero grade.

Examples
===
```
julia> x, y, z = basis((1,1,1))
3-element Vector{Blade{(1, 1, 1), Float64, UInt64, 1}}:
 1.0 v₁
 1.0 v₂
 1.0 v₃

julia> grade(x*y)
2

julia> grade(x + 1)
2-element Vector{Int64}:
 0
 1
```
"""
grade(::Type{<:HomogeneousMultivector{k}}) where k = k
grade(::HomogeneousMultivector{k}) where k = k
grade(a::MixedMultivector) = sort(unique(ublade_grade(u) for (u, v) ∈ comps(a) if !iszero(v)))

"""
`grade(a::AbstractMultivector, k)`

Grade projection of multivector `a` onto grade `k`. Returns a grade-`k` `Multivector`.

"""
grade(a::Blade, k) = grade(a) == k ? a : zero(a)
grade(a::Multivector{sig,C}, k) where {sig,C} = grade(a) == k ? a : zero(Multivector{sig,C,k})
function grade(a::MixedMultivector{sig,C}, k::Integer) where {sig,C}
	b = zero(Multivector{sig,C,k})
	for (u, v) ∈ comps(a)
		if ublade_grade(u) == k
			setcomp!(b, u, v)
		end
	end
	b
end

scalar(a::AbstractMultivector) = grade(a, 0)
scalar(a::Number) = a

isscalar(a::Blade) = iszero(grade(a))
isscalar(a::Multivector) = iszero(grade(a))
isscalar(a::MixedMultivector) = all(iszero.(grade(a)))

Base.isone(a::AbstractMultivector) = isscalar(a) && isone(a[])

# INDEXING

#=

Want to be able to index by integer or symbol (for labelled signatures),
and to support `OffsetSignature`s, for example, 0-based Lorentzian signatures.

=#

Base.getindex(a::AbstractMultivector) = getcomp(a, ublade_scalar(signature(a), keytype(a)))

#TODO: ambiguous, ugly, stupid
indices2ublade(::Any) = Int[]
indices2ublade(::Any, I::Integer...) = collect(I)
indices2ublade(::NamedTuple{labels}, I::Symbol...) where {labels} = [findfirst(==(i), labels) for i ∈ I]
indices2ublade(::Union{MetricSignature{sig},OffsetSignature{sig}}, I::Symbol...) where sig = indices2ublade(sig, I...) 
indices2ublade(::OffsetSignature{sig,indices}, I::Integer...) where {sig,indices} = [findfirst(==(i), indices) for i ∈ I]

function Base.getindex(a::AbstractMultivector, I...)
	ublade = indices2ublade(signature(a), I...)
	factor, ublade = ubladeprod(signature(a), ublade)
	ublade = convert_ublade(keytype(a), ublade)
	factor\getcomp(a, ublade)
end

function Base.setindex!(a::AbstractMultivector, v, I...)
	factor, ublade = ubladeprod(signature(a), collect(I))
	ublade = convert_ublade(keytype(a), ublade)
	setcomp!(a, ublade, factor\v)
end

# function Base.getindex(a::AbstractMultivector, I::Blade, Is::Blade...)
# 	blade = prod([I, Is...])
# 	blade.coeff\getcomp(a, blade.ublade)
# end