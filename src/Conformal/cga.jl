"""
# Conformal geometric algebra

Tools for working with conformal geometric algebra (CGA).

The algebra `CGA{Sig}` is an extension of the base space `Sig` with two additional dimensions
of square ``+1`` and ``-1``. Points in the base space `p::Multivector{Sig,1}` are associated to
null vectors in the higher space by
```julia
up(x) = n0 + x + x^2/2*noo
```
where `n0, noo = nullbasis(Sig)` are the null vectors representing the origin and the point at infinity.

Geometric primitives including _flats_ (points, point pairs, lines, planes, ...) and _rounds_
(point pairs, circles, spheres, ...) can be naturally represented as the inner or outer product
null spaces of CGA blades (see [`ipns`](@ref) and [`opns`](@ref)).

Geometric primitives may be translated in space with [`translate`](@ref) and combined with
the wedge product to produce meets (for `opns` geometries) or joins (intersections, for `ipns`).

Any CGA blade has a [`standardform`](@ref) from which geometric properties (such as position, radius
and carrier) can be easily retrieved.
"""
module Conformal

using StyledStrings
using ..GeometricAlgebra
using Compat

export CGA, CGABlade, CGAGeometry
export origin, infinity, nullbasis
export up, dn
export translate
export standardform
export ipns, opns

@compat public DirectionBlade, FlatBlade, DualFlatBlade, RoundBlade
@compat public FlatGeometry, RoundGeometry, PointAtInfinity, EmptySet


#= signature and basis display style =#

"""
Metric signature for the conformal geometric algebra over a base space with metric signature ``Sig``.
A conformal algebra `CGA{Sig}` has dimension `dimension(Sig) + 2`, with the two extra basis vectors squaring
to ``+1`` and ``-1``, respectively.
"""
abstract type CGA{Sig} end

GeometricAlgebra.dimension(::Type{CGA{Sig}}) where Sig = dimension(Sig) + 2
function GeometricAlgebra.basis_vector_square(P::Type{CGA{Sig}}, i::Integer) where Sig
	(GeometricAlgebra.canonical_signature(Sig)..., +1, -1)[i]
end
function GeometricAlgebra.get_basis_display_style(::Type{CGA{Sig}}) where Sig
	n = dimension(Sig)
	BasisDisplayStyle(n + 2, indices=[string.(1:n); 'p'; 'm'])
end


GeometricAlgebra.canonical_signature(::Type{CGA}) = (+1, -1)

#= signature promotion =#

GeometricAlgebra.signature_promote_rule(::Val{CGA{Sig}}, ::Val{Sig}) where Sig = CGA{Sig}
GeometricAlgebra.signature_convert(::Val{CGA{Sig}}, a::AbstractMultivector{Sig}) where Sig = embed(CGA{Sig}, a)


#= standard conformal null basis =#

nullbasis(S::Type{CGA{Sig}}) where Sig = (origin = origin(S), infinity = infinity(S))
origin(::Type{CGA{Sig}}) where Sig = Multivector{CGA{Sig},1}([zeros(dimension(Sig)); -0.5; 0.5])
infinity(::Type{CGA{Sig}}) where Sig = Multivector{CGA{Sig},1}([zeros(dimension(Sig)); +1; +1])

origin(n::Integer) = origin(CGA{n})
infinity(n::Integer) = infinity(CGA{n})
nullbasis(n::Integer) = nullbasis(CGA{n})

"""
	nullbasis(n) = (origin(n), infinity(n))
	origin(n)
	infinity(n)

Standard null basis vectors in the conformal geometric algebra `CGA{n}`.

The point at the origin ``e₀`` and the point at infinity ``e∞`` in the `n`-dimensional
conformal geometric algebra model are defined as
``e₀ = (e₊ + e₋)/2`` and ``e∞ = e₋ - e₊``
where ``e₊`` and ``e₋`` are the standard extra basis vectors squaring to ``+1`` and ``-1``.

Type-stable methods exist which accept the type `CGA{n}` instead of an integer `n`.
"""
nullbasis, origin, infinity



#= up/embedding maps =#

unembed(x::AbstractMultivector{CGA{Sig}}) where {Sig} = GeometricAlgebra.embed(Sig, x)

function up(x::Multivector{Sig,1}) where Sig
	o, oo = nullbasis(CGA{Sig})
	o + x + 2\(x⊙x)*oo
end
up(x::BasisBlade) = up(Multivector(x))
up(comps::AbstractVector) = up(Multivector{length(comps),1}(comps))
up(c, comps...) = up(Multivector{length(comps) + 1,1}(c, comps...))


function dn(x::Multivector{CGA{Sig},1}) where Sig
	oo = infinity(CGA{Sig})
	embed(Val(Sig), -scalar_prod(oo, x)\x)
end

"""
	up(::Multivector{Sig,1})::Multivector{CGA{Sig},1}
	dn(::Multivector{CGA{Sig},1})::Multivector{Sig,1}

"Lift up" a 1-vector in a base space `Sig` to a null vector in the 2d-up conformal algebra `CGA{Sig}`,
or "project down" a conformal 1-vector back into the base space.

The `up` map is given by
```math
up(x) = n0 + embed(x) + 1/2 x^2 noo
```
where `n0, noo = nullbasis(CGA{Sig})` are the points representing the origin and infinity.

For any vector `u` we have `dn(up(u)) == n` and `up∘dn` is idempotent.
"""
up, dn


#= translation operator =#

"""
	translate(p)
	translate(p, X) -> X′ = sandwich_prod(translate(p), X)

Translate `X::Multivector{CGA{Sig}}` by the displacement vector `p::Multivector{Sig,1}`.

The single-argument method returns the translation rotor and
the two-argument form applies the rotor to `X` with [`sandwich_prod`](@ref).

The translation rotor is defined as ``𝚃ₚ = \\exp(½ n_∞ p)`` where ``n_∞`` is
the point at [`infinity`](@ref).

# Examples
```jldoctest; setup = :(using GeometricAlgebra.Conformal)
julia> (p, x), noo = randn(Multivector{3,1}, 2), infinity(3);

julia> translate(p, up(x)) ≈ up(x + p)
true

julia> translate(p, noo) ≈ noo
true
```
"""
function translate(p::Grade{1,CGA{Sig}}) where Sig
	oo = infinity(CGA{Sig})
	1 + 2\oo∧p
end

translate(p::Grade{1,Sig}) where Sig = translate(embed(CGA{Sig}, p))
translate(p, X) = sandwich_prod(translate(p), X)



#= blade classification =#

"""
	CGABlade{Sig,K}:
		DirectionBlade{Sig,K}(E)
		FlatBlade{Sig,K}(E, p)
		DualFlatBlade{Sig,K}(E, p)
		RoundBlade{Sig,K}(E, p, r2)

Standard forms of blades in the conformal geometric algebra `CGA{Sig}` over base space `Sig`.

| Value | Mathematical form |
|:-----|:-----|
| `DirectionBlade(E)` | ``E ∧ n_∞`` |
| `FlatBlade(E, p)` | ``𝚃ₚ[n₀ ∧ E ∧ n_∞]`` |
| `DualFlatBlade(E, p)` | ``𝚃ₚ[E]`` |
| `RoundBlade(E, p, r2)` | ``𝚃ₚ[(n₀ + r2/2 n_∞) ∧ E]`` |

Any blade in `CGA{Sig}` is of exactly one of the forms above, where:
- ``E`` is a `K`-blade in the base space
- ``p`` is a position vector in the base space
- ``r2`` is a square-radius, which may be positive or negative
- ``n₀`` and ``n_∞`` are the points at the [`origin`](@ref) and at [`infinity`](@ref)
- ``Tₚ`` is the translation operator sending ``n₀`` to ``p``

The method [`standardform`](@ref) puts any blade in `CGA{Sig}` into one of these forms.
A standard blade `X::CGABlade` may be converted back to the usual additive form with `Multivector(X)`.

See table 14.1 of [^1] for discussion.

[^1]: Dorst, L., Fontijne, D., & Mann, S. (2010). Geometric Algebra for Computer Science: An Object-Oriented Approach to Geometry. Elsevier.
"""
abstract type CGABlade{Sig,K} end

struct DirectionBlade{Sig,K} <: CGABlade{Sig,K}
	E::Multivector{Sig,K}
end
struct FlatBlade{Sig,K} <: CGABlade{Sig,K}
	E::Multivector{Sig,K}
	p::Multivector{Sig,1}
end
struct DualFlatBlade{Sig,K} <: CGABlade{Sig,K}
	E::Multivector{Sig,K}
	p::Multivector{Sig,1}
end
struct RoundBlade{Sig,K} <: CGABlade{Sig,K}
	E::Multivector{Sig,K}
	p::Multivector{Sig,1}
	r2::Float64
end

@doc (@doc CGABlade) (DirectionBlade, FlatBlade, DualFlatBlade, RoundBlade)


GeometricAlgebra.Multivector(X::DirectionBlade{Sig}) where Sig = X.E ∧ infinity(Sig)
GeometricAlgebra.Multivector(X::FlatBlade{Sig}) where Sig = translate(X.p, origin(Sig) ∧ X.E ∧ infinity(Sig))
GeometricAlgebra.Multivector(X::DualFlatBlade{Sig}) where Sig = translate(X.p, X.E)
GeometricAlgebra.Multivector(X::RoundBlade{Sig}) where Sig = translate(X.p, (origin(Sig) + 2\X.r2*infinity(Sig)) ∧ X.E)

"""
	standardform(X::AbstractMultivector{CGA{Sig}}) -> CGABlade{Sig}

Put the blade `X` in standard form, returning a `CGABlade` object.
"""
function standardform(X::AbstractMultivector{<:CGA})
	@assert isblade(X)

	o = origin(signature(X))
	oo = infinity(signature(X))

	iszeroish(X) = isapprox(X, 0, atol=sqrt(eps(float(eltype(X)))))

	if iszeroish(X ∧ oo)
		if X ⨽ oo ≈ 0
			E = -unembed(X⨽o)
			DirectionBlade(E)
		else
			Xo = X⨽o # equal to -(o + p)∧E
			p = dn(X⨽Xo) # project origin onto X
			E = unembed(oo⨼Xo)
			FlatBlade(E, p)
		end
	else
		if iszeroish(X ⨽ oo)
			E = unembed(-(X ∧ oo)⨽o)
			if isscalar(E)
				DualFlatBlade(E, dn(o))
			else
				p = -E⨽unembed(X⨽o)/(E⊙E)
				DualFlatBlade(E, p)
			end
		else
			p = dn(sandwich_prod(X, oo))
			E = -involution(unembed(X ⨽ oo))
			r2 = (X⊙involution(X))/(E⊙E)
			RoundBlade(E, p, r2)
		end
	end
end
standardform(X::AbstractMultivector{Sig}) where Sig = standardform(embed(CGA{Sig}, X))


#= inner and outer product null spaces =#

abstract type CGAGeometry{Sig} end
struct FlatGeometry{Sig,K} <: CGAGeometry{Sig}
	p::Multivector{Sig,1}
	E::Multivector{Sig,K}
end
struct RoundGeometry{Sig,K} <: CGAGeometry{Sig}
	p::Multivector{Sig,1}
	E::Multivector{Sig,K}
	r2::Float64
end
struct PointAtInfinity{Sig} <: CGAGeometry{Sig} end
struct EmptySet{Sig} <: CGAGeometry{Sig} end

"""
	CGAGeometry{Sig}:
		FlatGeometry{Sig,K}(p, E)
		RoundGeometry{Sig,K}(p, E, r2)
		PointAtInfinity{Sig}()
		EmptySet{Sig}()

Subsets of ``ℝⁿ ∪ {∞}`` which are the [`ipns`](@ref) or [`opns`](@ref) of a conformal blade.

See also [`FlatGeometry`](@ref) and [`RoundGeometry`](@ref).
"""
CGAGeometry, PointAtInfinity, EmptySet

"""
	FlatGeometry{Sig,K}(p, E) <: CGAGeometry{Sig}

A `K`-flat in the base space `Sig` through the point `p::Multivector{Sig,1}`
spanning the `K`-blade `E::Multivector{Sig,K}`.

A `0`-flat is a point, a `1`-flat is a line, a `2`-flat is a plane, etc.
All flats include the unique point at infinity.

See also [`RoundGeometry`](@ref) and [`CGAGeometry`](@ref).
"""
FlatGeometry

"""
	RoundGeometry{Sig,K}(p, E, r2) <: CGAGeometry{Sig}

A `K`-round in the base space `Sig` around the point `p::Multivector{Sig,1}`
spanning the `K`-blade `E::Multivector{Sig,K}` with square radius `r2`.

A `0`-round is a point, a `1`-round is a point pair, a `2`-round is a circle,
a `3`-round is a sphere, etc. Rounds never include the point at infinity.

!!! note
	The square radius `r2` may be negative, in which case the round is formally
	the empty set, but one may also interpret this as an "imaginary" radius.

See also [`FlatGeometry`](@ref) and [`CGAGeometry`](@ref).
"""
RoundGeometry

ipns(X::DirectionBlade{Sig}) where Sig = PointAtInfinity{Sig}()
ipns(X::DualFlatBlade) = FlatGeometry(X.p, hodgedual(X.E))
ipns(X::FlatBlade{Sig}) where Sig = EmptySet{Sig}()
ipns(X::RoundBlade) = RoundGeometry(X.p, hodgedual(X.E), -X.r2)

opns(X::DirectionBlade{Sig}) where Sig = PointAtInfinity{Sig}()
opns(X::DualFlatBlade{Sig}) where Sig = EmptySet{Sig}()
opns(X::FlatBlade) = FlatGeometry(X.p, X.E)
opns(X::RoundBlade) = RoundGeometry(X.p, X.E, X.r2)

ipns(X::AbstractMultivector) = ipns(standardform(X))
opns(X::AbstractMultivector) = opns(standardform(X))

"""
	ipns(A) -> CGAGeometry
	opns(A) -> CGAGeometry

Inner or outer product null space of the blade `A`.

This is the set of points ``x`` satisfying `up(x)⋅A ≈ 0` (IPNS) or `up(x)∧A ≈ 0` (OPNS),
possibly including the point at infinity.
"""
ipns, opns



#= display methods =#

function showfields(io::IO, X::T) where T
	iszero(nfields(X)) && return
	pad = maximum(length.(string.(fieldnames(T))))
	for field in fieldnames(T)
		printstyled(io, "\n  ", rpad(field, pad), color=:cyan)
		print(io, " = ")
		val = getfield(X, field)
		if val isa Multivector
			GeometricAlgebra.show_multivector(io, val, inline=true, showzeros=false)
		else
			show(io, val)
		end
	end
end

GeometricAlgebra.grade(X::DirectionBlade) = grade(X.E) + 1
GeometricAlgebra.grade(X::FlatBlade) = grade(X.E) + 2
GeometricAlgebra.grade(X::DualFlatBlade) = grade(X.E)
GeometricAlgebra.grade(X::RoundBlade) = grade(X.E) + 1

showformula(::Type{<:DirectionBlade}) = "E∧oo"
showformula(::Type{<:FlatBlade}) = "translate(p, n0∧E∧noo)"
showformula(::Type{<:DualFlatBlade}) = "translate(p, E)"
showformula(::Type{<:RoundBlade}) = "translate(p, (n0 + r2/2*noo)∧E)"
function Base.show(io::IO, mime::MIME"text/plain", X::T) where T <: CGABlade{Sig,K} where {Sig,K}
	println(io, T, ":")
	print(io, " $(grade(X))-blade of the form ")
	printstyled(io, showformula(T), color=:cyan)
	print(io, ":")
	showfields(io, X)
end

showformula(::Type{<:PointAtInfinity}) = ""
showformula(::Type{<:EmptySet}) = ""
function showformula(::Type{<:FlatGeometry{Sig,K}}) where {Sig,K}
	desc = get(("point", "line", "plane", "volume"), K + 1, nothing)
	desc = isnothing(desc) ? "" : " ($desc)"
	styled"$K-flat$desc through {cyan:p} spanning {cyan:E}"
end
function showformula(::Type{<:RoundGeometry{Sig,K}}) where {Sig,K}
	desc = get(("point", "point pair", "circle", "sphere"), K + 1, nothing)
	desc = isnothing(desc) ? "" : " ($desc)"
	styled"$K-round$desc around center {cyan:p} spanning {cyan:E} with square radius {cyan:r2}"
end
function Base.show(io::IO, mime::MIME"text/plain", X::T) where T <: Union{FlatGeometry,RoundGeometry}
	println(io, T, ": ")
	print(io, " ", showformula(T), ":")
	showfields(io, X)
end




end # module Conformal