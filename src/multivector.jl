#= MULTIVECTORS

There are three types representing objects in the geometric algebra, all
subtypes of `AbstractMultivector`.
 - `Blade`: a scalar multiple of a wedge product of basis vectors
 - `Multivector`: a homogeneous multivector; a sum of same-grade blades
 - `MixedMultivector`: an inhomogeneous multivector. All objects in the algebra
	can be expressed as a mixed multivector.
=#

"""
	AbstractMultivector{sig,C}

Supertype of all elements in the geometric algebra over a vector space with
metric signature `sig`.

The parameter `C` is the type which the multivector components are internally stored as.
For instance, `a::AbstractMultivector{(1,1), SparseVector{Int64}}` is a (mixed)
multivector in the Euclidean plane with components of `eltype(a) == Int64`
stored in a `Vector`.
"""
abstract type AbstractMultivector{sig,C} <: Number end


"""
	Blade{sig,k,T,B} <: AbstractMultivector{sig,Pair{B,T}}
	Blade{sig}(coeff, ublade)

A grade `k` blade (i.e., wedge product of `k` basis vectors) with coefficient of type `T`
and unit blade represented by type `B`, in vector space of metric signature `sig`.

Unit blade type `B` | E.g. for ``v_1∧v_3∧v_4`` | Signature
:-------------------|:------------------|:----------------
Unsigned            | `0b1101`          | any
Vector{<:Integer}   | `[1, 3, 4]`       | any
Vector{<:Symbol}    | `[:v1, :v2, :v3]` | labelled

Examples
===
```jldoctest
julia> Blade{(1,1,1)}(42, 0b101)
Grade-2 Blade{⟨+++⟩, 2, Int64, UInt8}:
 42 v13

julia> Blade{(x=1,y=1,z=1)}(1, [:x])
Grade-1 Blade{⟨x+,y+,z+⟩, 1, Int64, Array{Symbol,1}}:
 1 x
```
"""
struct Blade{sig,k,T,B} <: AbstractMultivector{sig,Pair{B,T}}
	coeff::T
	ublade::B
	function Blade{sig,k,T,B}(coeff, ublade) where {sig,k,T,B}
		ublade = convert_ublade(sig, B, ublade)
		# @assert ublade_grade(ublade) == k
		new{sig,k,T,B}(coeff, ublade)
	end
end

Blade{sig}(coeff::T, ublade::B) where {sig,T,B} = Blade{sig,ublade_grade(ublade),T,B}(coeff, ublade)
Blade{sig,k}(coeff::T, ublade::B) where {sig,k,T,B} = Blade{sig,k,T,B}(coeff::T, ublade::B)

function (::Type{Blade{sig,k,T,B} where k})(coeff, ublade=ublade_scalar(B)) where {sig,T,B}
	ublade = convert(B, ublade)
	k = ublade_grade(ublade)
	Blade{sig,k,T,B}(coeff, ublade)
end



"""
	Multivector{sig,k,C} <: AbstractMultivector{sig,C}

Homogeneous multivector of grade `k` in space of metric signature `sig`, containing
`binomial(dimension(sig), k)` components.

Components are stored in the field `.comps` of type `C`, which in general may be an
`AbstractVector{T}` or an `AbstractDict{B,T}`, where `T` is the component type.

If `C<:AbstractVector`, components are stored in lexicographic order.
E.g., the component of a 3-multivector corresponding to ``v_1∧v_3∧v_4``
or `0b1101` is located at the `ublade2lindex(0b1101) == 3`rd index.
"""
struct Multivector{sig,k,C} <: AbstractMultivector{sig,C}
	comps::C
	function Multivector{sig,k,C}(comps) where {sig,k,C<:AbstractVector}
		# @assert length(comps) == binomial(dimension(sig), k)
		new{sig,k,C}(comps)
	end
	function Multivector{sig,k,C}(comps) where {sig,k,C<:AbstractDict}
		# @assert all(ublade_grade(u) == k for u ∈ keys(comps))
		new{sig,k,C}(comps)
	end
end
Multivector{sig,k}(comps::C) where {sig,k,C} = Multivector{sig,k,C}(comps)



"""
	MixedMultivector{sig,C} <: AbstractMultivector{sig,C}
	MixedMultivector{sig}(comps)

Inhomogeneous multivector in vector space of metric signature `sig`,
containing `2^dimension(sig)` components of varying grade.

Components are stored in the field `.comps` of type `C`, which in
general may be an `AbstractVector{T}` or an `AbstractDict{B,T}`,
where `T` is the component type (similar to `Multivectors{sig,k,C}`).

When `C<:AbstractVector`, components are ordered by the binary value of the unit blade.
E.g., the coefficient of ``v_1∧v_3∧v_4`` is stored at index `1 + Int(0b1101) == 14`.

Examples
===
```jldoctest
julia> MixedMultivector{(1,1,1)}(Dict(Int[] => 1, [1] => 2, [1,2] => 3))
MixedMultivector{⟨+++⟩, Dict{Array{Int64,1},Int64}}:
 1
 2 v1
 3 v12

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
	function MixedMultivector{sig,C}(comps::AbstractVector) where {sig,C<:AbstractVector}
		@assert length(comps) == 2^dimension(sig)
		new{sig,C}(comps)
	end
	function MixedMultivector{sig,C}(comps::AbstractDict) where {sig,C<:AbstractDict}
		new{sig,C}(comps)
	end
end
MixedMultivector{sig}(comps::C) where {sig,C} = MixedMultivector{sig,C}(comps)


# USEFUL TYPE ALIASES

const Scalar = Union{filter(T -> !(T <: AbstractMultivector), subtypes(Number))...}
const HomogeneousMultivector{k} = Union{Blade{sig,k},Multivector{sig,k}} where sig
const CompositeMultivector{C} = Union{Multivector{sig,k,C},MixedMultivector{sig,C}} where {sig,k}
const multivector_types = [Blade, Multivector, MixedMultivector]




"""
	eltype(::AbstractMultivector)
	eltype(::Type{<:AbstractMultivector})

Give the numerical type of the components of a multivector.
"""
Base.eltype

Base.eltype(::Type{<:Blade{sig,k,T} where {sig,k}}) where T = T
Base.eltype(::Type{<:CompositeMultivector{C}}) where {C} = valtype(C)

Base.valtype(a::Union{Type{<:AbstractMultivector},<:AbstractMultivector}) = eltype(a)

"""
	keytype(::AbstractMultivector)
	keytype(::Type{<:AbstractMultivector})

Give the type in which unit blades are encoded in a multivector.

E.g., for `Blade{sig,k,T,B}`, this is `B`.
For multivectors with components stored in a dictionary of type `C`, this is `keytype(C)`.
For components stored in a vector (where the unit blades are not stored explicitly),
`keytype(a)` defaults to an unsigned type used to index the vector.

Examples
===
```jldoctest
julia> Blade{(1,1,1)}(42, 0b101)
Grade-2 Blade{⟨+++⟩, 2, Int64, UInt8}:
 42 v13

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


narrowest_uint(n) = n <= 8 ? UInt8 : n <= 16 ? UInt16 : n <= 32 ? UInt32 : n <= 64 ? UInt64 : n <= 128 ? UInt128 : Vector{Int}
best_ublade_type(sig) = sig_has_dim(sig) ? narrowest_uint(dimension(sig)) : Vector{Int}
Base.keytype(::Type{<:Blade{sig,k,T,B} where k}) where {sig,T,B} = B
Base.keytype(::Type{<:CompositeMultivector{<:AbstractVector}}) = UInt # appropriate?
Base.keytype(::Type{<:CompositeMultivector{<:AbstractDict{B}}}) where B = B

Base.keytype(::T) where T<:AbstractMultivector = keytype(T)



"""
	signature(::AbstractMultivector)
	signature(::Type{<:AbstractMultivector})

Return the metric signature associated with the given multivector or multivector type. 

Metric signatures can be any object implementing `getindex` to return the norm of a
given basis vector. Examples include `Tuple`s, `NamedTuple`s, or instances of `AbstractMetricSignature`.

Examples
===
```jldoctest
julia> x = Blade{(x=1, y=1)}(1, [:x])
Grade-1 Blade{⟨x+,y+⟩, 1, Int64, Array{Symbol,1}}:
 1 x

julia> signature(x)
(x = 1, y = 1)
```
"""
signature(::Type{<:AbstractMultivector{sig}}) where sig = sig
signature(::AbstractMultivector{sig}) where sig = sig

"""
	dimension(::AbstractMultivector)
	dimension(::Type{<:AbstractMultivector})

Return the dimension of the vector space underlying the given multivector
(or multivector type).
"""
dimension(::Type{<:AbstractMultivector{sig}}) where sig = dimension(sig)
dimension(::AbstractMultivector{sig}) where sig = dimension(sig)



# BLADE ITERATION

"""
	blades(a::AbstractMultivector)

Return a vector of `Blades` representing the non-zero components of the multivector `a`.

Examples
===

```jldoctest
julia> x, y, z = basis((1,1,1));

julia> blades(1 + 2x + 3y*z)
3-element Array{Blade{(1, 1, 1),k,Float64,UInt64} where k,1}:
 1.0
 2.0 v1
 3.0 v23
```
"""
blades(a::Blade) = iszero(a) ? typeof(a)[] : [a]
function blades(a::Multivector{sig,k,<:AbstractVector}) where {sig,k}
	T = Blade{sig,k,eltype(a),keytype(a)} # for type stability
	T[Blade{sig}(coeff, ublade) for (coeff, ublade) ∈ zip(a.comps, FixedGradeBlades{keytype(a)}(k, dimension(sig))) if !iszero(coeff)]
end
function blades(a::MixedMultivector{sig,<:AbstractVector}) where {sig}
	T = Blade{sig,k,eltype(a),keytype(a)} where k # for type stability
	T[T(coeff, convert_ublade(a, unsigned(ublade - 1))) for (ublade, coeff) ∈ enumerate(a.comps) if !iszero(coeff)]
end
function blades(a::Multivector{sig,k,<:AbstractDict}) where {sig,k}
	T = Blade{sig,k,eltype(a),keytype(a)} # for type stability
	T[T(coeff, ublade) for (ublade, coeff) ∈ a.comps if !iszero(coeff)]
end
function blades(a::MixedMultivector{sig,<:AbstractDict}) where {sig}
	T = Blade{sig,k,eltype(a),keytype(a)} where k # for type stability
	T[T(coeff, ublade) for (ublade, coeff) ∈ a.comps if !iszero(coeff)]
end




# COMPONENT ACCESS

getcomp(a::Blade{sig,k,T,B}, ublade::B) where {sig,k,T,B} = a.ublade == ublade ? a.coeff : zero(T)
function getcomp(a::Multivector{sig,k,<:AbstractVector}, ublade::Unsigned) where {sig,k}
	ublade_grade(ublade) == k || return zero(eltype(a))
	i = ublade2lindex(ublade)
	a.comps[i]
end
getcomp(a::MixedMultivector{sig,<:AbstractVector}, ublade::Unsigned) where sig = a.comps[begin + ublade]
getcomp(a::CompositeMultivector{<:AbstractDict{B}}, ublade::B) where {B} = get(a.comps, ublade, zero(eltype(a)))
getcomp(a::AbstractMultivector, ublade) = getcomp(a, convert_ublade(a, ublade))


setcomp!(::Blade, ublade, v) = error("cannot set component of blade")
function setcomp!(a::Multivector{sig,k,<:AbstractVector}, ublade, v) where {sig,k}
	ublade_grade(ublade) == k || error("cannot set non-$k grade component of $k-multivector")
	isempty(a.comps) && error("$k-multivector in $(dimension(a)) dimensions has zero components")
	i = ublade2lindex(ublade)
	a.comps[i] = v
end
setcomp!(a::MixedMultivector{sig,<:AbstractVector}, ublade::Unsigned, v) where sig = a.comps[begin + ublade] = v
setcomp!(a::CompositeMultivector{<:AbstractDict{B}}, ublade::B, v) where {B} = a.comps[ublade] = v
setcomp!(a::CompositeMultivector, ublade, v) = setcomp!(a, convert_ublade(a, ublade), v)


function add!(a::AbstractMultivector, b::Scalar)
	ublade = ublade_scalar(keytype(a))
	setcomp!(a, ublade, getcomp(a, ublade))
	a
end
function add!(a::AbstractMultivector, b::AbstractMultivector)
	# iszero(b) && return a # if b is empty, blades(b) is empty too
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
	Blade{signature(a)}(f(a), a.ublade)
end
mapcomps(f::Function, a::AbstractMultivector) = mapcomps!(f, deepcopy(a))



## ADDITIVE/MULTIPLICATIVE IDENTITIES

_zeros(::Type{Vector{T}}, n) where T = zeros(T, n)
_zeros(::Type{SparseVector{T}}, n) where T = spzeros(T, n)

Base.zero(::Type{<:Blade{sig,k,T,B}}) where {sig,k,T,B} = Blade{sig,k,T,B}(zero(T), ublade_first_of_grade(B, k))
Base.zero(::Type{<:Multivector{sig,k,C}}) where {sig,k,C<:AbstractVector} = Multivector{sig,k}(_zeros(C, binomial(dimension(sig), k)))
Base.zero(::Type{<:MixedMultivector{sig,C}}) where {sig,C<:AbstractVector} = MixedMultivector{sig}(_zeros(C, 2^dimension(sig)))
Base.zero(T::Type{<:CompositeMultivector{C}}) where {C<:AbstractDict} = T(C())

# for non-concrete types
Base.zero(::Type{<:Blade{sig,k,T}}) where {sig,k,T} = zero(Blade{sig,k,T,UInt})
Base.zero(::Type{<:Blade{sig,k}}) where {sig,k} = zero(Blade{sig,k,Float64})
Base.zero(::Type{<:Multivector{sig,k}}) where {sig,k} = zero(Multivector{sig,k,Vector{Float64}})
Base.zero(::Type{<:MixedMultivector{sig}}) where {sig} = zero(MixedMultivector{sig,Vector{Float64}})

Base.iszero(a::Blade) = iszero(a.coeff)
Base.iszero(a::CompositeMultivector{<:AbstractVector}) = iszero(a.comps)
Base.iszero(a::CompositeMultivector{<:AbstractDict}) = all(iszero(v) for (u, v) ∈ a.comps)


Base.one(::Type{<:Blade{sig,k,T,B}}) where {sig,k,T,B} = Blade{sig}(one(T), ublade_scalar(B))
Base.oneunit(a::Blade{sig,k,T,B}) where {sig,k,T,B} = Blade{sig}(one(T), a.ublade)

function Base.one(T::Type{<:Multivector{sig,k,C}}) where {sig,k,C}
	a = zero(Multivector{sig,0,C}) # scalar / 0-grade
	a[] = one(eltype(a))
	a
end
function Base.one(T::Type{<:MixedMultivector})
	a = zero(T)
	a[] = one(eltype(a))
	a
end

Base.zero(a::AbstractMultivector) = zero(typeof(a)) 
Base.one(a::AbstractMultivector) = one(typeof(a))





#= MULTIVECTOR TYPE GYMNASTICS

Different `<:AbstractMultivector{sig}` instances are _compatible_ if they
share the same metric signature `sig`. Compatible multivectors may be promoted
to a common type.
This involves promoting the `eltype`, ublade type, and for `CompositeMultivector`s,
the storage type (whether components are stored in a vector, dict, etc).
Promotion does *not* change the grade parameter `k` of `HomogeneousMultivector`s.
=#

# assert that multivectors share same metric signature and return it
shared_sig(::Type{<:AbstractMultivector{sig}}...) where sig = sig
shared_sig(::Type{<:AbstractMultivector}...) = error("multivectors must share the same metric signature")
shared_sig(as::AbstractMultivector...) = shared_sig(typeof.(as)...)



# set the grade parameter `k` of a (homogeneous) multivector type, dispatching on `k::Integer`
set_grade(::Type{<:Blade{sig,k′,T,B} where k′}, ::Val{k}) where {sig,T,B,k} = Blade{sig,k,T,B}
set_grade(::Type{<:Multivector{sig,k′,C} where k′}, ::Val{k}) where {sig,C,k} = Multivector{sig,k,C}
set_grade(::Type{<:Blade{sig,k′,T,B} where k′}, ::Missing) where {sig,T,B} = Blade{sig,k,T,B} where k
set_grade(::Type{<:Multivector{sig,k′,C} where k′}, ::Missing) where {sig,C} = Multivector{sig,k,C} where k
set_grade(T::Type{<:MixedMultivector}, ::Missing) = T

set_grade(T, k::Integer) = set_grade(T, Val(k)) # not type stable

set_eltype(::Type{<:Vector}, T) = Vector{T}
set_eltype(::Type{<:Dict{B}}, T) where B = Dict{B,T}

# note: should work for UnionAll's with parametric k
set_eltype(::Type{<:Blade{sig,k,T′,B} where T′}, T) where {sig,k,B}   = Blade{sig,k,T,B}
set_eltype(::Type{<:Blade{sig,k,T′,B} where {k,T′}}, T) where {sig,B} = Blade{sig,k,T,B} where k
set_eltype(::Type{<:Multivector{sig,k,C}}, T) where {sig,k,C}         = Multivector{sig,k,set_eltype(C, T)}
set_eltype(::Type{<:Multivector{sig,k,C} where k}, T) where {sig,C}   = Multivector{sig,k,set_eltype(C, T)} where k
set_eltype(::Type{<:MixedMultivector{sig,C}}, T) where {sig,C}        = MixedMultivector{sig,set_eltype(C, T)}

set_storagetype(::Type{<:Multivector{sig,k}}, C) where {sig,k} = Multivector{sig,k,C}
set_storagetype(::Type{<:MixedMultivector{sig}}, C) where {sig} = MixedMultivector{sig,C}

# `StorageType` is a wrapper for multivector storage types (e.g., `Vector{Float64}` or `Dict{<ublade_type>,<eltype>}`)
# which harnesses the built-in promotion system.
struct StorageType{T}
	# StorageType(T) = new{T}()
end
unwrap(::Type{StorageType{T}}) where T = T

# for `Blade`s, `storagetype` returns the preferred default storage type for composite multivectors.
storagetype(::Type{<:Blade{sig,k,T,B} where k}) where {sig,T,B<:Union{UInt8,UInt16,UInt32,UInt64}} = StorageType{Vector{T}}
storagetype(::Type{<:Blade{sig,k,T,B} where k}) where {sig,T,B} = StorageType{Dict{B,T}}
storagetype(::Type{<:CompositeMultivector{C}}) where {C} = StorageType{C}

# as storage types, dicts are more general than vectors
Base.promote_rule(::Type{StorageType{Vector{T}}}, ::Type{StorageType{Vector{S}}}) where {T,S} = StorageType{Vector{promote_type(T, S)}}
Base.promote_rule(::Type{StorageType{Dict{B,T}}}, ::Type{StorageType{Vector{S}}}) where {B,T,S} = StorageType{Dict{B,promote_type(T, S)}}
Base.promote_rule(::Type{StorageType{Dict{B1,T1}}}, ::Type{StorageType{Dict{B2,T2}}}) where {B1,B2,T1,T2} = StorageType{Dict{promote_ublade_type(B1, B2),promote_type(T1, T2)}}

	
"""
	_best_type(::Type{M}, as::Type{<:AbstractMultivector}...)

Given `M` in `(Blade, Multivector, MixedMultivector)`, give the best concrete subtype of `M`
which can simultaneously represent objects with types `as...`.
E.g., `_best_type(M, Blade{:sig,2,UInt8,Int8}, Multivector{:sig,1,Vector{Int},Float64})`
returns a concrete `M` type with eltype `Float64` and ublade type `Vector{Int}`.
"""
function _best_type(::Type{Blade}, as::Type{<:AbstractMultivector}...)
	sig = shared_sig(as...)
	T = promote_type(eltype.(as)...)
	B = promote_type(keytype.(as)...)
	Blade{sig,k,T,B} where k
end
function _best_type(::Type{Multivector}, as::Type{<:AbstractMultivector}...)
	sig = shared_sig(as...)
	C = promote_type(storagetype.(as)...) |> unwrap
	Multivector{sig,k,C} where k
end
function _best_type(::Type{MixedMultivector}, as::Type{<:AbstractMultivector}...)
	sig = shared_sig(as...)
	C = promote_type(storagetype.(as)...) |> unwrap
	MixedMultivector{sig,C}
end

# called with a single argument, returns the same type
# and preserves grade (unlike e.g., _best_type(Blade, ::Blade), which leaves `k` free)
_best_type(T::Type{<:AbstractMultivector}) = set_grade(T, missing)

# @generated to achieve type stability
@generated _best_type(::Type{M}, a::AbstractMultivector, bs::AbstractMultivector...) where M = _best_type(M, a, bs...)


function best_type(args...; grade=missing, el=nothing)
	T = _best_type(args...)
	isnothing(el) || (T = set_eltype(T, promote_type(eltype(T), el)))
	set_grade(T, grade)
end


#= CONSTRUCTORS TO CONVERT FROM LOWER MULTIVECTOR TYPES

As sets, Blade ⊂ Multivector ⊂ MixedMultivector.
Conversion into more general types is possible.
=#

Blade(a::Blade) = a
Multivector(a::HomogeneousMultivector{k}) where k = add!(zero(best_type(Multivector, a; grade=Val(k))), a)
MixedMultivector(a::AbstractMultivector) = add!(zero(best_type(MixedMultivector, a)), a)



## CONVERSION
# need to be able to convert to types returned by `best_type`
# must convert to type identically except for possibly the grade parameter (which should not be changed)

# Conversion of element type
Base.convert(::Type{<:Blade{sig,k,T,B} where k}, a::Blade) where {sig,T,B} = Blade{sig,grade(a),T,B}(a.coeff, convert_ublade(sig, B, a.ublade))
Base.convert(::Type{<:Multivector{sig,k,C} where k}, a::Multivector) where {sig,C} = add!(zero(Multivector{sig,grade(a),C}), a)
Base.convert(::Type{<:MixedMultivector{sig,C}}, a::MixedMultivector) where {sig,C} = add!(zero(MixedMultivector{sig,C}), a)

# Conversion from lower multivector type
Base.convert(T::Type{<:Multivector{sig,k,C}}, a::Blade) where {sig,k,C} = add!(zero(Multivector{sig,k,C}), a)
Base.convert(T::Type{<:Multivector{sig,k,C} where k}, a::Blade) where {sig,C} = add!(zero(Multivector{sig,grade(a),C}), a)
Base.convert(T::Type{<:MixedMultivector}, a::HomogeneousMultivector) = add!(zero(T), a)

# Conversion from scalar
Base.convert(::Type{<:Blade{sig,k,T,B} where k}, a::Scalar) where {sig,T,B} = Blade{sig,0,T,B}(a, ublade_scalar(B))
Base.convert(::Type{<:Multivector{sig,k,C} where k}, a::Scalar) where {sig,C} = zero(Multivector{sig,0,C}) + a
Base.convert(T::Type{<:MixedMultivector}, a::Scalar) = zero(T) + a


## PROMOTION
# promotion should give same type except for the grade parameter, which may remain different.

# Same multivector type
Base.promote_rule(T::Type{<:Blade}, S::Type{<:Blade}) = best_type(Blade, T, S)
Base.promote_rule(T::Type{<:Multivector}, S::Type{<:Multivector}) = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector}, S::Type{<:MixedMultivector}) = best_type(MixedMultivector, T, S)

# Promotion to higher multivector type
Base.promote_rule(T::Type{<:Multivector}, S::Type{<:Blade}) = best_type(Multivector, T, S)
Base.promote_rule(T::Type{<:MixedMultivector}, S::Type{<:AbstractMultivector}) = best_type(MixedMultivector, T, S)

# Promotion from scalars
Base.promote_rule(T::Type{<:AbstractMultivector}, S::Type{<:Scalar}) = best_type(T; el=S)

# hack: need to treat `<:HomogeneousMultivector`s of nearly-identical types except for
# differing grade as "the same type" as far as promotion is concerned
Base.promote(as::(Blade{sig,k,T,B} where k)...) where {sig,T,B} = as
Base.promote(as::(Multivector{sig,k,C} where k)...) where {sig,C} = as


grade_maybe(::Type{<:HomogeneousMultivector{k}}) where k = Val(k)
grade_maybe(::Type{<:AbstractMultivector}) = missing

# set_eltype(T::Type{<:HomogeneousMultivector{k}}, S) where k = best_type(T; el=S, grade=Val(k))
# set_eltype(T::Type{<:MixedMultivector}, S) = best_type(T; el=Val(S))
# set_eltype(T::Type{<:AbstractMultivector}, S) = best_type(T; el=S, grade=grade_maybe(T))

Base.float(T::Type{<:AbstractMultivector}) = set_eltype(T, float(eltype(T)))
Base.float(a::AbstractMultivector) = convert(float(typeof(a)), a)

Base.big(T::Type{<:AbstractMultivector}) = set_eltype(T, big(eltype(T)))
Base.big(a::AbstractMultivector) = convert(big(typeof(a)), a)

Base.complex(T::Type{<:AbstractMultivector}) = set_eltype(T, complex(eltype(T)))
Base.complex(a::AbstractMultivector) = convert(complex(typeof(a)), a)

_constructor(::Multivector{sig,k}) where {sig,k} = Multivector{sig,k}
_constructor(::MixedMultivector{sig}) where {sig} = MixedMultivector{sig}

Base.real(a::Blade{sig,k,T,B}) where {sig,k,T,B} = Blade{sig,k,real(T),B}(real(a.coeff), a.ublade)
Base.real(a::CompositeMultivector{<:AbstractVector}) = _constructor(a)(real(a.comps))
Base.real(a::CompositeMultivector{<:AbstractDict}) = _constructor(a)(mapcomps(a, real))

Base.imag(a::Blade{sig,k,T,B}) where {sig,k,T,B} = Blade{sig,k,real(T),B}(imag(a.coeff), a.ublade)
Base.imag(a::CompositeMultivector{<:AbstractVector}) = _constructor(a)(imag(a.comps))
Base.imag(a::CompositeMultivector{<:AbstractDict}) = _constructor(a)(mapcomps(a, imag))



## EQUALITY

==(a::T, b::T) where T<:Blade = a.ublade == b.ublade && a.coeff == b.coeff
==(a::T, b::T) where T<:Multivector = a.comps == b.comps
==(a::T, b::T) where T<:MixedMultivector = a.comps == b.comps

# ==(::AbstractMultivector{sig1}, ::AbstractMultivector{sig2}) where {sig1,sig2} = false

==(a::Blade{sig,k1,T,B}, b::Blade{sig,k2,T,B}) where {sig,k1,k2,T,B} = iszero(a) && iszero(b)
==(a::Multivector{sig,k1,C}, b::Multivector{sig,k2,C}) where {sig,k1,k2,C} = iszero(a) && iszero(b)
==(::MixedMultivector{sig1,C}, ::MixedMultivector{sig2,C}) where {sig1,sig2,C} = false

==(a::AbstractMultivector, b::Scalar) = isscalar(a) && scalar(a) == b
==(a::Scalar, b::AbstractMultivector) = isscalar(b) && scalar(b) == a


## APPROXIMATE EQUALITY

function _isapprox(a::(Blade{sig,k,T,B} where k), b::(Blade{sig,k,T,B} where k); kwargs...) where {sig,T,B}
	if a.ublade == b.ublade
		isapprox(a.coeff, b.coeff; kwargs...)
	else # orthogonal blades can be approximately equal if they are both approximately zero
		isapprox(a.coeff, zero(T); kwargs...) && isapprox(b.coeff, zero(T); kwargs...)
	end
end

# between composite multivectors of identical type (and grade)
function _isapprox(a::T, b::T; kwargs...) where T<:CompositeMultivector{<:AbstractVector}
	isapprox(a.comps, b.comps; kwargs...)
end
function _isapprox(a::T, b::T; kwargs...) where T<:CompositeMultivector{<:AbstractDict}
	B = keytype(T)
	for i ∈ 1:2^dimension(T)
		ublade = convert_ublade(sig, B, unsigned(i - 1))
		isapprox(getcomp(a, ublade), getcomp(b, ublade); kwargs...) || return false
	end
	true
end

function Base.isapprox(a::HomogeneousMultivector, b::HomogeneousMultivector; kwargs...)
	if grade(a) == grade(b)
		_isapprox(promote(a, b)...; kwargs...)
	else # objects of differing grade can still be approximately equal if they both approximately vanish
		@show typeof(zero(a))
		_isapprox(a, zero(a); kwargs...) && _isapprox(b, zero(b); kwargs...)
	end
end

Base.isapprox(a::AbstractMultivector, b::AbstractMultivector; kwargs...) = _isapprox(promote(a, b)...; kwargs...)
Base.isapprox(a::AbstractMultivector, b::Scalar; kwargs...) = _isapprox(promote(a, b)...; kwargs...)
Base.isapprox(a::Scalar, b::AbstractMultivector; kwargs...) = _isapprox(promote(a, b)...; kwargs...)






"""
	vol(x)

The volume form or pseudoscalar element.

Examples
===
```jldoctest
julia> x = basis((1, 1, 1), 1)
Grade-1 Blade{⟨+++⟩, 1, Float64, UInt64}:
 1.0 v1

julia> vol(x)
Grade-3 Blade{⟨+++⟩, 3, Float64, UInt64}:
 1.0 v123
```
"""
function vol(T::Type{<:AbstractMultivector{sig}}) where sig
	Blade{sig}(one(eltype(T)), ublade_first_of_grade(keytype(T), dimension(sig)))
end
vol(::T) where T<:AbstractMultivector = vol(T)




# GRADE SELECTION

"""
	grade(a::AbstractMultivector)

The grade of the homogeneous multivector `a`.
For `MixedMultivector`s, see `grades`.

Examples
===
```jldoctest
julia> x, y, z = basis((1,1,1))
3-element Array{Blade{⟨+++⟩, 1, Float64, UInt64},1}:
 1.0 v1
 1.0 v2
 1.0 v3

julia> grade(x*y)
2
```
"""
grade(x::Scalar, k) = iszero(k) ? x : zero(x)
grade(::Type{<:HomogeneousMultivector{k}}) where k = k
grade(::HomogeneousMultivector{k}) where k = k
grade(a::MixedMultivector) = error("mixed multivectors do not have a grade; use grades(::MixedMultivector) instead")

"""
	grades(a::AbstractMultivector)

The non-zero grades of the multivector `a`.

Examples
===
```jldoctest
julia> x, y, z = basis((1,1,1))
3-element Array{Blade{⟨+++⟩, 1, Float64, UInt64},1}:
 1.0 v1
 1.0 v2
 1.0 v3

julia> grades(1 + x*y)
2-element Array{Int64,1}:
 0
 2
```
"""
grades(::HomogeneousMultivector{k}) where k = [k]
grades(a::MixedMultivector) = sort(unique(grade(u) for u ∈ blades(a) if !iszero(u)))

"""
	grade(a::AbstractMultivector, k)

The grade `k` part of a (mixed) multivector `a`. Returns a grade `k` `Multivector`.

Examples
===
```jldoctest; setup = :( (x, y) = basis((x=1,y=1)) )
julia> grade(x + 1 + y + x*y, 1)
Grade-1 Multivector{⟨x+,y+⟩, 1, Array{Float64,1}}:
 1.0 x
 1.0 y
```
"""
grade(a::Blade{sig,j}, k) where {sig,j} = j == k ? a : zero(a)
grade(a::Multivector{sig,j}, k) where {sig,j} = j == k ? a : zero(a)
function grade(a::MixedMultivector{sig,C}, k::Integer) where {sig,C}
	b = zero(Multivector{sig,k,C})
	for u ∈ blades(a)
		if grade(u) == k
			add!(b, u)
		end
	end
	b
end

scalar(a::Number) = a
scalar(a::AbstractMultivector) = grade(a, 0)[]

# Good idea?
# pseudoscalar(a::Number) = zero(a)
# pseudoscalar(a::AbstractMultivector) = first(blades(grade(a, dimension(a)))).coeff

isscalar(a::Number) = true
isscalar(a::HomogeneousMultivector{0}) = true
isscalar(a::HomogeneousMultivector{k}) where k = iszero(k) ? true : iszero(a)
isscalar(a::MixedMultivector) = all(iszero.(grades(a)))

Base.isone(a::AbstractMultivector) = isscalar(a) && isone(scalar(a))

# INDEXING

#=

Want to be able to index by integer or symbol (for labelled signatures),
and to support `OffsetSignature`s, for example, 0-based Lorentzian signatures.

=#

Base.getindex(a::AbstractMultivector) = getcomp(a, ublade_scalar(keytype(a)))
Base.setindex!(a::AbstractMultivector, v) = setcomp!(a, ublade_scalar(keytype(a)), v)

# convert e.g., [:t, :y] or [0, 2] to canonical [1, 3]
function interpret_ublade(a::AbstractMultivector{sig}, I) where sig
	ublade = [unoffset_index(sig, i) for i ∈ I] # not necessarily canonical
	factor, ublade = ubladeprod(sig, ublade)
	ublade = convert_ublade(a, ublade)
	factor, ublade
end

function _getindex(a::AbstractMultivector, I...)	
	factor, ublade = interpret_ublade(a, I)
	factor\getcomp(a, ublade)
end

function _setindex!(a::AbstractMultivector, v, I::Union{<:Integer,Symbol}...)
	factor, ublade = interpret_ublade(a, I)
	setcomp!(a, ublade, factor\v)
end

Base.getindex(a::AbstractMultivector, I::Integer...) = _getindex(a, I...)
Base.getindex(a::AbstractMultivector, I::Symbol...) = _getindex(a, I...)
Base.setindex!(a::AbstractMultivector, v, I::Integer...) = _setindex!(a, v, I...)
Base.setindex!(a::AbstractMultivector, v, I::Symbol...) = _setindex!(a, v, I...)

function Base.getindex(a::AbstractMultivector, I::Blade...)
	blade = prod(collect(I))
	blade.coeff\getcomp(a, blade.ublade)
end
function Base.setindex!(a::AbstractMultivector, v, I::Blade...)
	blade = prod(collect(I))
	setcomp!(a, blade.ublade, blade.coeff\v)
end
