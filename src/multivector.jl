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

Supertype of all elements in the geometric algebra over a vector space with
metric signature `sig`.

The parameter `C` is the type which the multivector components are stored as.
For instance, `a::AbstractMultivector{(1,1), SparseVector{Int64}}` is a (mixed)
multivector in the Euclidean plane with components of `eltype(a) == Int64`
stored in a `Vector`.
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

A grade `k` blade (i.e., wedge product of `k` basis vectors) with coefficient of type `T`
and unit blade represented by type `B`, in vector space of metric signature `sig`.

Unit blade type `B` | E.g. for ``v_1∧v_3∧v_4`` | Signature
:-------------------|:------------------|:----------------
Unsigned            | `0b1101`          | any
Vector{<:Integer}   | `[1, 3, 4]`       | any
Vector{<:Symbol}    | `[:v₁, :v₃, :v₄]`  | labelled

Examples
===
```jldoctest
julia> Blade{(1,1,1)}(42, 0b101)
2-grade Blade{⟨+++⟩, Int64, UInt8, 2}:
 42 v₁v₃

julia> Blade{(x=1,y=1,z=1)}(1, [:x])
1-grade Blade{⟨x+,y+,z+⟩, Int64, Array{Symbol,1}, 1}:
 1 x
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

Components are stored in the field `.comps` of type `C`, which in general may be an
`AbstractVector{T}` or an `AbstractDict{B,T}`, where `T` is the component type.

If `C<:AbstractVector`, components are stored in lexicographic order.
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

Multivector{sig,C}(args...) where {sig,C} = error("must initialize multivector with grade parameter")

"""
```
Multivector{sig}(k, comps)
Multivector{sig}(comps::AbstractVector)
```

Construct a multivector of grade `k` from a vector or dict. When the argument `k`
is omitted, returns a 1-multivector (if `length(comps) == dim(sig)`).

Examples
===
```jldoctest
julia> Multivector{(1,1,1)}([10, 0, 20])
1-grade Multivector{⟨+++⟩, Array{Int64,1}, 1}:
 10 v₁
 20 v₃

julia> Multivector{(x=1,y=1)}(2, [pi])
2-grade Multivector{⟨x+,y+⟩, Array{Irrational{:π},1}, 2}:
 π xy
```
"""
Multivector{sig}(comps::C) where {sig,C<:AbstractVector} = Multivector{sig,C,1}(comps)

"""
```
MixedMultivector{sig,C} <: AbstractMultivector{sig,C}
MixedMultivector{sig}(comps)
```

Inhomogeneous multivector in vector space of metric signature `sig`,
containing `2^dim(sig)` independent components of mixed grade.

Components are stored in the field `.comps` of type `C`, which in
general may be an `AbstractVector{T}` or an `AbstractDict{B,T}`,
where `T` is the component type (similar to `Multivectors{sig,C,k}`).

When `C<:AbstractVector`, components are ordered by the binary value of the unit blade.
E.g., the coefficient of `v1∧v3∧v4` is stored at index `1 + Int(0b1101) == 14`.

Examples
===
```jldoctest
julia> MixedMultivector{(1,1,1)}(Dict(Int[] => 1, [1] => 2, [1,2] => 3))
MixedMultivector{⟨+++⟩, Dict{Array{Int64,1},Int64}}:
 1
 2 v₁
 3 v₁v₂

julia> basis((x=1, y=1), 1) + 7
MixedMultivector{⟨x+,y+⟩, Array{Float64,1}}:
 7.0
 1.0 x

julia> ans.comps
4-element Array{Float64,1}:
 7.0
 1.0
 0.0
 0.0
```

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





"""
	eltype(a::AbstractMultivector)

Give the type of the components of the multivector `a`.
"""
Base.eltype

"""
	keytype(a::AbstractMultivector)

Give the type in which unit blades are encoded in the multivector `a`.

For an instance of `Blade{sig,T,B,k}`, this is `B`.
For multivectors with components stored in a dictionary of type `C`, this is `keytype(C)`.
For components stored in a vector (where the unit blades are not stored explicitly),
`keytype(a)` defaults to an unsigned type.

Examples
===
```jldoctest
julia> Blade{(1,1,1)}(42, 0b101)
2-grade Blade{⟨+++⟩, Int64, UInt8, 2}:
 42 v₁v₃

julia> keytype(ans)
UInt8

julia> MixedMultivector{EuclideanSignature}(Dict([:z] => 1, [:x, :y] => 2))
MixedMultivector{EuclideanSignature, Dict{Array{Symbol,1},Int64}}:
 1 z
 2 xy

julia> keytype(ans)
Array{Symbol,1}
```
"""
Base.keytype


Base.eltype(::Type{<:Blade{sig,T}}) where {sig,T} = T
Base.eltype(::Type{<:CompositeMultivector{sig,<:AbstractVector{T}}}) where {sig,T} = T
Base.eltype(::Type{<:CompositeMultivector{sig,<:AbstractDict{B,T}}}) where {sig,B,T} = T

narrowest_uint(n) = n <= 8 ? UInt8 : n <= 16 ? UInt16 : n <= 32 ? UInt32 : n <= 64 ? UInt64 : n <= 128 ? UInt128 : Vector{Int}
best_ublade_type(sig) = sig_has_dimension(sig) ? narrowest_uint(dim(sig)) : Vector{Int}
Base.keytype(::Type{<:Blade{sig,T,B}}) where {sig,T,B} = B
Base.keytype(::Type{<:CompositeMultivector{sig,<:AbstractVector}}) where {sig} = best_ublade_type(sig)
Base.keytype(::Type{<:CompositeMultivector{sig,<:AbstractDict{B}}}) where {sig,B} = B
Base.keytype(a::AbstractMultivector) = keytype(typeof(a))




"""
`signature(::AbstractMultivector)`
`signature(::Type{<:AbstractMultivector})`

Return the metric signature associated with the given multivector or multivector type. 

Metric signatures can be any object implementing `getindex` to return the norm of a
given basis vector.
Signatures may be `Tuple`s, `NamedTuple`s, or instances of `AbstractMetricSignature`.
"""
signature(::Type{<:AbstractMultivector{sig}}) where sig = sig
signature(a::AbstractMultivector) = signature(typeof(a))

"""
	dim(a)

Return the dimension of the vector space underlying the given multivector
(or multivector type).
"""
dim(a::Union{Type{AbstractMultivector},AbstractMultivector}) = dim(signature(a))





# COMPONENT ACCESS

getcomp(a::Blade{sig,T,B}, ublade::B) where {sig,T,B} = a.ublade == ublade ? a.coeff : zero(T)
function getcomp(a::Multivector{sig,<:AbstractVector,k}, ublade::Unsigned) where {sig,k}
	ublade_grade(ublade) == k || return zero(eltype(a))
	i = ublade2lindex(ublade)
	a.comps[i]
end
getcomp(a::MixedMultivector{sig,<:AbstractVector}, ublade::Unsigned) where sig = a.comps[begin + ublade]
getcomp(a::CompositeMultivector{sig,<:AbstractDict{B}}, ublade::B) where {sig,B} = get(a.comps, ublade, zero(eltype(a)))
getcomp(a::AbstractMultivector, ublade) = getcomp(a, convert_ublade(a, ublade))


setcomp!(::Blade, ublade, v) = error("cannot set component of blade")
function setcomp!(a::Multivector{sig,<:AbstractVector,k}, ublade, v) where {sig,k}
	ublade_grade(ublade) == k || error("cannot set non-$k grade components of $k-multivector")
	i = ublade2lindex(ublade)
	a.comps[i] = v
end
setcomp!(a::MixedMultivector{sig,<:AbstractVector}, ublade::Unsigned, v) where sig = a.comps[begin + ublade] = v
setcomp!(a::CompositeMultivector{sig,<:AbstractDict{B}}, ublade::B, v) where {sig,B} = a.comps[ublade] = v
setcomp!(a::CompositeMultivector, ublade, v) = setcomp!(a, convert_ublade(a, ublade), v)



"""
	blades(a::AbstractMultivector)

Return a vector of `Blades` representing the non-zero components of the multivector `a`.

Examples
===

```jldoctest
julia> x, y, z = basis((1,1,1));


julia> blades(1 + 2x + 3y*z)
3-element Array{Blade{(1, 1, 1),Float64,UInt8,k} where k,1}:
 1.0
 2.0 v₁
 3.0 v₂v₃
```
"""
blades(a::Blade) = iszero(a) ? typeof(a)[] : [a]
blades(a::Multivector{sig,<:AbstractVector,k}) where {sig,k} =
	[Blade{sig}(coeff, ublade) for (coeff, ublade) ∈ zip(a.comps, FixedGradeBlades{keytype(a)}(k, dim(sig))) if !iszero(coeff)]
blades(a::MixedMultivector{sig,<:AbstractVector}) where {sig} =
	[Blade{sig}(coeff, convert_ublade(a, unsigned(ublade - 1))) for (ublade, coeff) ∈ enumerate(a.comps) if !iszero(coeff)]
blades(a::CompositeMultivector{sig,<:AbstractDict}) where sig =
	[Blade{sig}(coeff, ublade) for (ublade, coeff) ∈ a.comps if !iszero(coeff)]

blade_bvs(a::Blade{sig,T,Vector{<:Integer}}) where {sig,T} = a.ublade
blade_bvs(a::Blade{sig,T,<:Unsigned}) where {sig,T} = convert_ublade(sig, Vector{Int}, a.ublade)



function add!(a::AbstractMultivector, b::AbstractMultivector)
	iszero(b) && return a
	for u ∈ blades(b)
		setcomp!(a, u.ublade, getcomp(a, u.ublade) + u.coeff)
	end
	a
end


"""
	mapcomps!(f, a)
	mapcomps(f, a)

Map the components of a multivector `a` by the function `f`, which
is called with a single `Blade` argument and whose return type should be `eltype(a)`.

The mutating version `mapcomps!` maps the components in-place.
"""
mapcomps!, mapcomps

function mapcomps!(f::Function, a::CompositeMultivector)
	for u ∈ blades(a)
		setcomp!(a, u.ublade, f(u))
	end
	a
end
function mapcomps(f::Function, a::Blade)
	Blade{signature(a)}(f(a.ublade, a.coeff), a.ublade)
end
mapcomps(f::Function, a::AbstractMultivector) = mapcomps!(f, deepcopy(a))



## ADDITIVE/MULTIPLICATIVE IDENTITIES

_zeros(::Type{Vector{T}}, n) where T = zeros(T, n)
_zeros(::Type{SparseVector{T}}, n) where T = spzeros(T, n)

Base.zero(::Type{<:Blade{sig,T,B,k}}) where {sig,T,B,k} = Blade{sig,T,B,k}(zero(T), ublade_first_of_grade(B, k))
Base.zero(::Type{<:Multivector{sig,C,k}}) where {sig,C<:AbstractVector,k} = Multivector{sig}(k, _zeros(C, binomial(dim(sig), k)))
Base.zero(::Type{<:MixedMultivector{sig,C}}) where {sig,C<:AbstractVector} = MixedMultivector{sig}(_zeros(C, 2^dim(sig)))
Base.zero(T::Type{<:CompositeMultivector{sig,C}}) where {sig,C<:AbstractDict} = T(C())

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector{sig,C}) where {sig,C<:AbstractVector} = iszero(a.comps)
Base.iszero(a::CompositeMultivector{sig,C}) where {sig,C<:AbstractDict} = all(iszero(v) for (u, v) ∈ a.comps)


Base.one(::Type{<:Blade{sig,T,B}}) where {sig,T,B} = Blade{sig}(one(T), ublade_scalar(B))
Base.oneunit(a::Blade{sig,T,B,k}) where {sig,T,B,k} = Blade{sig}(one(T), a.ublade)

function Base.one(T::Type{<:Multivector{sig,C}}) where {sig,C}
	a = zero(Multivector{sig,C,0}) # scalar / 0-grade
	setcomp!(a, ublade_scalar(keytype(T)), one(eltype(T)))
	a
end
function Base.one(T::Type{<:MixedMultivector})
	a = zero(T)
	setcomp!(a, ublade_scalar(keytype(T)), one(eltype(T)))
	a
end

# not needed if AbstractMultivector <: Number
Base.zero(a::AbstractMultivector) = zero(typeof(a)) 
Base.one(a::AbstractMultivector) = one(typeof(a))





# MULTIVECTOR TYPE GYMNASTICS

# assert that multivectors share same metric signature and return it
shared_sig(::Type{<:AbstractMultivector{sig}}...) where sig = sig
shared_sig(::Type{<:AbstractMultivector}...) = error("multivectors must share the same metric signature")
shared_sig(as::AbstractMultivector...) = shared_sig(typeof.(as)...)

# give the most sensible storage type of a CompositeMultivector with given signature, eltype, and keytype
storagetype(sig, T, B::Type{<:Unsigned}) = dim(sig) <= 8 ? Vector{T} : SparseVector{T}
storagetype(sig, T, B) = Dict{B,T}


# give the most sensible AbstractMultivector type which can represent the given arguments
function best_type(M::Type{Blade}, as::Type{<:AbstractMultivector}...; k=nothing, T=Union{})
	sig = shared_sig(as...)
	S = promote_type(T, eltype.(as)...)
	B = promote_ublade_type(keytype.(as)...)
	isnothing(k) ? Blade{sig,S,B} : Blade{sig,S,B,k}
end
function best_type(M::Type{Multivector}, as::Type{<:AbstractMultivector}...; k=nothing, T=Union{})
	sig = shared_sig(as...)
	S = promote_type(T, eltype.(as)...)
	B = promote_ublade_type(keytype.(as)...)
	C = storagetype(sig, S, B)
	isnothing(k) ? Multivector{sig,C} : Multivector{sig,C,k}
end
function best_type(M::Type{MixedMultivector}, as::Type{<:AbstractMultivector}...; T=Union{})
	sig = shared_sig(as...)
	S = promote_type(T, eltype.(as)...)
	B = promote_ublade_type(keytype.(as)...)
	C = storagetype(sig, S, B)
	MixedMultivector{sig,C}
end
best_type(M, as::AbstractMultivector...; args...) = best_type(M, typeof.(as)...; args...)



#= CONVERSION

As sets, Blade ⊂ Multivector ⊂ MixedMultivector.
Conversion into more general types is possible.
=#

Base.convert(::Type{T}, a::T) where T<:AbstractMultivector = a

## conversions preserving storage type

function Base.convert(::Type{Multivector}, a::HomogeneousMultivector{k}) where k
	b = zero(best_type(Multivector, a; k))
	for u ∈ blades(a)
		add!(b, u)
	end
	b
end

function Base.convert(::Type{MixedMultivector}, a::AbstractMultivector)
	b = zero(best_type(MixedMultivector, a))
	for u ∈ blades(a)
		add!(b, u)
	end
	b
end

## conversions with specified target signature and storage type

function Base.convert(::Type{<:Blade{sig,T,B}}, a::Blade) where {sig,T,B}
	Blade{sig}(convert(T, a.coeff), convert_ublade(sig, B, a.ublade))
end

function Base.convert(::Type{<:Multivector{sig,C}}, a::HomogeneousMultivector{k}) where {sig,C,k}
	b = zero(Multivector{sig,C,k})
	for u ∈ blades(a)
		add!(b, u)
	end
	b
end

function Base.convert(::Type{MixedMultivector{sig,C}}, a::AbstractMultivector) where {sig,C}
	b = zero(MixedMultivector{sig,C})
	for u ∈ blades(a)
		add!(b, u)
	end
	b
end


## PROMOTION

Base.promote_rule(T::Type{<:Blade}, S::Type{<:Blade}) = best_type(Blade, T, S)
Base.promote_rule(T::Type{<:Multivector}, S::Type{<:HomogeneousMultivector}) = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector}, S::Type{<:AbstractMultivector}) = best_type(MixedMultivector, T, S)





## EQUALITY

for eq ∈ [:(==), :(Base.isapprox)]
	#TODO: should 1e-18 v₁ == 1e-18 v₂ ?
	@eval $eq(a::T, b::T; kwargs...) where T<:Blade = a.ublade == b.ublade && $eq(a.coeff, b.coeff; kwargs...)
	@eval $eq(a::Multivector{sig,<:AbstractVector,k}, b::Multivector{sig,<:AbstractVector,k}; kwargs...) where {sig,k} = $eq(a.comps, b.comps; kwargs...)
	@eval $eq(a::MixedMultivector{sig,<:AbstractVector}, b::MixedMultivector{sig,<:AbstractVector}; kwargs...) where sig = $eq(a.comps, b.comps; kwargs...)

	@eval function $eq(a::T, b::T; kwargs...) where T<:CompositeMultivector{sig,<:AbstractDict} where sig
			@assert keytype(a) == keytype(b)
			for u ∈ union(keys(a.comps), keys(b.comps))
				$eq(getcomp(a, u), getcomp(b, u); kwargs...) ? continue : return false
			end
			true
		end

	# fallback: coerce into same type
	@eval $eq(a::AbstractMultivector, b::AbstractMultivector; kwargs...) = $eq(promote(a, b)...; kwargs...)

	@eval $eq(a::AbstractMultivector, b::Number) = isscalar(a) && $eq(a[], b)
	@eval $eq(a::Number, b::AbstractMultivector) = b == a
end






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
basis(M::Type{Blade{sig,T,B,k}}, i::Integer) where {sig,T,B,k} = M(one(T), convert_ublade(M, lindex2ublade(keytype(M), k, i)))
function basis(M::Type{Multivector{sig,C,k}}, i::Integer) where {sig,C<:AbstractVector,k}
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
basis(sig, i::Integer) = basis(Blade{sig,Float64,best_ublade_type(sig),1}, i)
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
vol(T::Type{<:AbstractMultivector{sig}}) where sig = Blade{sig}(one(eltype(T)), ublade_first_of_grade(keytype(T), dim(sig)))
vol(::T) where T<:AbstractMultivector = vol(T)




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
grade(a::MixedMultivector) = sort(unique(grade(u) for u ∈ blades(a) if !iszero(u)))

"""
`grade(a::AbstractMultivector, k)`

Grade projection of multivector `a` onto grade `k`. Returns a grade-`k` `Multivector`.

"""
grade(a::Blade, k) = grade(a) == k ? a : zero(best_type(Blade, a; k=0))
grade(a::Multivector{sig,C}, k) where {sig,C} = grade(a) == k ? a : zero(Multivector{sig,C,k})
function grade(a::MixedMultivector{sig,C}, k::Integer) where {sig,C}
	b = zero(Multivector{sig,C,k})
	for u ∈ blades(a)
		if grade(u) == k
			add!(b, u)
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

Base.getindex(a::AbstractMultivector) = getcomp(a, ublade_scalar(keytype(a)))
Base.setindex!(a::AbstractMultivector, v) = setcomp!(a, ublade_scalar(keytype(a)), v)

function Base.getindex(a::AbstractMultivector, I...)
	ublade = collect(I)
	factor, ublade = ubladeprod(signature(a), collect(I))
	ublade = convert_ublade(a, ublade)
	factor\getcomp(a, ublade)
end

function Base.setindex!(a::AbstractMultivector, v, I...)
	factor, ublade = ubladeprod(signature(a), collect(I))
	ublade = convert_ublade(a, ublade)
	setcomp!(a, ublade, factor\v)
end

# function Base.getindex(a::AbstractMultivector, I::Blade, Is::Blade...)
# 	blade = prod([I, Is...])
# 	blade.coeff\getcomp(a, blade.ublade)
# end