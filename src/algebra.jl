#= Equality =#

Base.:(==)(a::BasisBlade{Sig}, b::BasisBlade{Sig}) where Sig = bitsof(a) == bitsof(b) ? a.coeff == b.coeff : iszero(a) && iszero(b)
Base.:(==)(a::KVector{Sig}, b::KVector{Sig}) where {Sig} = grade(a) == grade(b) ? a.comps == b.comps : iszero(a) && iszero(b)
Base.:(==)(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig} = a.comps == b.comps

Base.:(==)(a::AbstractMultivector, b::Number) = isscalar(a) && scalar(a) == b
Base.:(==)(a::Number, b::AbstractMultivector) = isscalar(b) && a == scalar(b)

# equality between different multivector types
# TODO: implement without conversions?
Base.:(==)(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig} = let T = largest_type(a, b)
	T(a) == T(b)
end



#= Approximate Equality =#

isapproxzero(a; kwargs...) = isapprox(a, zero(a); kwargs...)
isapproxzero(a::BasisBlade; kwargs...) = isapproxzero(a.coeff; kwargs...)
isapproxzero(a::CompositeMultivector; kwargs...) = isapproxzero(a.comps; kwargs...)

Base.isapprox(a::BasisBlade{Sig}, b::BasisBlade{Sig}; kwargs...) where Sig = bitsof(a) == bitsof(b) ? isapprox(a.coeff, b.coeff; kwargs...) : isapproxzero(a) && isapproxzero(b)
Base.isapprox(a::KVector{Sig}, b::KVector{Sig}; kwargs...) where {Sig} = grade(a) == grade(b) ? isapprox(a.comps, b.comps; kwargs...) : isapproxzero(a) && isapproxzero(b)
Base.isapprox(a::Multivector{Sig}, b::Multivector{Sig}; kwargs...) where {Sig} = isapprox(a.comps, b.comps; kwargs...)

# promote scalar to target multivector type and compare component arrays
Base.:isapprox(a::BasisBlade, b::Number; kwargs...) = isapprox(KVector(a), b; kwargs...)
Base.:isapprox(a::CompositeMultivector, b::Number; kwargs...) = isapprox(a, zero(a) + b; kwargs...)
Base.:isapprox(a::Number, b::AbstractMultivector; kwargs...) = isapprox(b, a; kwargs...)

Base.isapprox(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}; kwargs...) where {Sig} = let T = largest_type(a, b)
	isapprox(T(a), T(b); kwargs...)
end



#= Scalar Multiplication =#

const Scalar = Union{Number,SymbolicUtils.Symbolic}

scalar_multiply(a::BasisBlade{Sig,K}, b) where {Sig,K} = BasisBlade{Sig,K}(bitsof(a) => a.coeff*b)
scalar_multiply(a, b::BasisBlade{Sig,K}) where {Sig,K} = BasisBlade{Sig,K}(bitsof(b) => a*b.coeff)

scalar_multiply(a::CompositeMultivector, b) = constructor(a)(a.comps*b)
scalar_multiply(a, b::CompositeMultivector) = constructor(b)(a*b.comps)

Base.:*(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
Base.:*(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)
Base.:-(a::AbstractMultivector) = -numberone(eltype(a))*a

promote_to(T, x) = convert(promote_type(T, typeof(x)), x)
Base.:/(a::AbstractMultivector, b::Scalar) = a*inv(promote_to(eltype(a), b))
Base.:\(a::Scalar, b::AbstractMultivector) = inv(promote_to(eltype(b), a))*b

Base.://(a::AbstractMultivector, b::Scalar) = a*(one(b)//b)



#= Addition =#

add!(a::KVector, b::KVector) = (a.comps .+= b.comps; a)
add!(a::Multivector, b::Multivector) = (a.comps .+= b.comps; a)

add!(a::CompositeMultivector, b::BasisBlade) = (a.comps[componentindex(a, b)] += b.coeff; a)
add!(a::CompositeMultivector, b::CompositeMultivector) = (a.comps[componentindex(a, b)] += b.comps; a)


# add alike types
Base.:+(As::KVector{Sig,K}...) where {Sig,K} = KVector{Sig,K}(sum(a.comps for a ∈ As))
Base.:+(As::Multivector{Sig}...) where {Sig} = Multivector{Sig}(sum(a.comps for a ∈ As))

# convert unalike to alike # TODO: reduce intermediate allocations
Base.:+(As::HomogeneousMultivector{Sig,K}...) where {Sig,K} = +(KVector.(As)...)
Base.:+(As::AbstractMultivector{Sig}...) where {Sig} = +(Multivector.(As)...)

Base.:-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)



#= Scalar Addition =#

add_scalar(a::AbstractMultivector{Sig}, b) where {Sig} = a + BasisBlade{Sig,0}(0 => b)
function add_scalar(a::Multivector, b)
	constructor(a)(copy_setindex(a.comps, a.comps[1] + b, 1))
end

Base.:+(a::AbstractMultivector, b::Scalar) = add_scalar(a, b)
Base.:+(a::Scalar, b::AbstractMultivector) = add_scalar(b, a)

Base.:-(a::AbstractMultivector, b::Scalar) = add_scalar(a, -b)
Base.:-(a::Scalar, b::AbstractMultivector) = add_scalar(-b, a)



#= Geometric Multiplication =#

# compute multivector type suitable to represent the result of f(a, b)
# leaves eltype/storage type parameters variable
result_type(f::Any, a::AbstractMultivector, b::AbstractMultivector) = result_type(f, typeof(a), typeof(b))

geometric_prod(a::Scalar, b::Scalar) = a*b
geometric_prod(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
geometric_prod(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)

function geometric_prod(a::BasisBlade{Sig}, b::BasisBlade{Sig}) where {Sig}
	factor, bits = geometric_prod_bits(Sig, bitsof(a), bitsof(b))
	BasisBlade{Sig}(bits => factor*(a.coeff*b.coeff))
end

# subspace of blades is closed under *
result_type(::typeof(geometric_prod), ::Type{<:BasisBlade{Sig}}, ::Type{<:BasisBlade{Sig}}) where {Sig} = BasisBlade{Sig}
function result_type(::typeof(geometric_prod), ::Type{<:HomogeneousMultivector{Sig,P}}, ::Type{<:HomogeneousMultivector{Sig,Q}}) where {Sig,P,Q}
	k, K = minmax(P, Q)
	k == 0 && return KVector{Sig,K} # multiplication by scalar preserves grade
	K == dimension(Sig) && return KVector{Sig,dimension(Sig) - k} # multiplication by pseudoscalar flips grade
	Multivector{Sig} # product of homogeneous multivectors is in general inhomogeneous
end
result_type(::typeof(geometric_prod), ::Type{<:AbstractMultivector{Sig}}, ::Type{<:AbstractMultivector{Sig}}) where {Sig} = Multivector{Sig}

function _geometric_prod(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	c = zero(similar(result_type(geometric_prod, a, b), a, b))
	for (abits::UInt, acoeff) ∈ nonzero_components(a), (bbits::UInt, bcoeff) ∈ nonzero_components(b)
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		i = componentindex(c, bits)
		c = setindex!!(c, c.comps[i] + factor*(acoeff*bcoeff), i)
	end
	c
end

@generated geometric_prod(a, b) = symbolic_optim(_geometric_prod, a, b)

Base.:*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)



#= Scalar Product =#

"""
	a ⊙ b
	scalar_prod(a, b) -> Number

Scalar part of the multivector product `a*b`.
"""
scalar_prod(a::Scalar, b::Scalar) = a*b
scalar_prod(a::AbstractMultivector, b::Scalar) = scalar(a)*b
scalar_prod(a::Scalar, b::AbstractMultivector) = a*scalar(b)

scalar_prod(a::BasisBlade{Sig,K}, b::BasisBlade{Sig,K}) where {Sig,K} = bitsof(a) == bitsof(b) ? scalar(a*b) : numberzero(promote_type(eltype(a), eltype(b)))
scalar_prod(a::BasisBlade{Sig}, b::BasisBlade{Sig}) where {Sig} = numberzero(promote_type(eltype(a), eltype(b)))

# these are 𝒪(n) while scalar(geom_prod(a, b)) is 𝒪(n^2)
scalar_prod(a::KVector{Sig,K}, b::KVector{Sig,K}) where {Sig,K} = sum(geometric_square_factor.(Ref(Sig), bitsof(a)) .* (a.comps .* b.comps))
scalar_prod(a::KVector{Sig}, b::KVector{Sig}) where {Sig} = numberzero(promote_type(eltype(a), eltype(b)))

scalar_prod(a::Multivector{Sig}, b::KVector{Sig,K}) where {Sig,K} = scalar_prod(grade(a, K), b)
scalar_prod(a::KVector{Sig,K}, b::Multivector{Sig}) where {Sig,K} = scalar_prod(a, grade(b, K))

scalar_prod(a::CompositeMultivector{Sig}, b::BasisBlade{Sig}) where {Sig} = geometric_square_factor(Sig, bitsof(b))*(a.comps[componentindex(a, b)]*b.coeff)
scalar_prod(a::BasisBlade{Sig}, b::CompositeMultivector{Sig}) where {Sig} = geometric_square_factor(Sig, bitsof(a))*(a.coeff*b.comps[componentindex(b, a)])

scalar_prod(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig} = sum(geometric_square_factor.(Ref(Sig), bitsof(a)) .* (a.comps .* b.comps))

"$(@doc scalar_prod)"
a ⊙ b = scalar_prod(a, b)



#= Graded Products =#

"""
	graded_prod(grade_selector::Function, a, b)

A "graded" product of multivectors, generalising the wedge ``∧``, inner ``⋅`` and contraction products.
For example, the wedge product is defined as:
```julia
wedge(a, b) = graded_prod(+, a, b)
```

If `grade(a) == p` and `grade(b) == q`, then `graded_prod(f, a, b)` is the grade `f(p, q)` part of `a*b`.
The definition extends linearly to general multivectors ``A`` and ``B`` as
```math
	(A, B) ↦ \\sum_{p,q} ⟨⟨A⟩_p ⟨B⟩_q⟩_{f(p, q)}
```
where ``⟨⋅⟩_k`` denotes the grade ``k`` part.
"""
graded_prod(grade_selector, a::Scalar, b::Scalar) = iszero(grade_selector(0, 0)) ? a*b : zero(promote_type(a, b))
graded_prod(grade_selector, a::AbstractMultivector, b::Scalar) = grade(a) == (grade_selector(grade(a), 0)) ? scalar_multiply(a, b) : zero(a)
graded_prod(grade_selector, a::Scalar, b::AbstractMultivector) = grade(b) == (grade_selector(0, grade(b))) ? scalar_multiply(a, b) : zero(a)

function graded_prod(grade_selector::Function, a::BasisBlade{Sig}, b::BasisBlade{Sig}) where {Sig}
	if count_ones(bitsof(a) ⊻ bitsof(b)) == grade_selector(grade(a), grade(b))
		a*b
	else
		BasisBlade{Sig}(0 => zero(promote_type(eltype(a), eltype(b))))
	end
end

result_type(::Tuple{typeof(graded_prod),Any}, ::Type{<:BasisBlade{Sig}}, ::Type{<:BasisBlade{Sig}}) where {Sig} = BasisBlade{Sig}
result_type(::Tuple{typeof(graded_prod),GradeSelector}, ::Type{<:HomogeneousMultivector{Sig,P}}, ::Type{<:HomogeneousMultivector{Sig,Q}}) where {Sig,P,Q,GradeSelector} = KVector{Sig,GradeSelector.instance(P, Q)}
result_type(::Tuple{typeof(graded_prod),GradeSelector}, ::Type{<:AbstractMultivector{Sig}}, ::Type{<:AbstractMultivector{Sig}}) where {Sig,GradeSelector} = Multivector{Sig}

function _graded_prod(grade_selector::Function, a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	c = zero(similar(result_type((graded_prod, grade_selector), a, b), a, b))
	for (abits, acoeff) ∈ nonzero_components(a), (bbits, bcoeff) ∈ nonzero_components(b)
		bits = abits ⊻ bbits
		if count_ones(bits) == grade_selector(count_ones(abits), count_ones(bbits))
			factor = geometric_prod_factor(Sig, abits, bbits)
			i = componentindex(c, bits)
			c = setindex!!(c, c.comps[i] + factor*(acoeff*bcoeff), i)
		end
	end
	c
end

@generated function graded_prod(grade_selector, a, b)
	symbolic_optim(a, b) do a, b
		_graded_prod(grade_selector.instance, a, b)
	end
end



#= Derived Products =#

"""
	a ∧ b
	wedge(a, b)

Wedge product of multivectors, a.k.a. the _outer_, _exterior_ or _alternating_ product.

This is a grade-raising operation, equivalent to [`graded_prod(+, a, b)`](@ref).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ∧ b` is the grade ``p + q`` part of `a*b`.
"""
wedge(a, b) = graded_prod(+, a, b)

"$(@doc wedge)"
a ∧ b = wedge(a, b)

"""
	a ⋅ b
	inner(a, b)

Inner product of multivectors.

This is a grade lowering operation, equivalent to [`graded_prod(abs∘-, a, b)`](@ref).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ⋅ b` is the grade ``|p - q|`` part of `a*b`.

Note that for scalars `a` and `b`, the inner product reduces to scalar multiplication,
in contrast to some authors (see [^D02] for discussion).
"""
inner(a, b) = graded_prod(abs∘-, a, b)

"""
$(@doc inner)
[^D02]: Leo Dorst, "The Inner Products of Geometric Algebra", 2002. [doi:10.1007/978-1-4612-0089-5_2](https://dx.doi.org/10.1007/978-1-4612-0089-5_2)
"""
⋅(a, b) = inner(a, b)

"""
	a ⨼ b
	lcontract(a, b)

Left contraction of multivectors. See also [`rcontract`](@ref).

Equivalent to [`graded_prod((p, q) -> q - p, a, b)`](@ref).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ⨼ b` is the grade ``q - p`` part of `a*b`.
"""
lcontract(a, b) = graded_prod((-)∘-, a, b)

"$(@doc lcontract)"
⨼(a, b) = lcontract(a, b)

"""
	a ⨽ b
	rcontract(a, b)

Left contraction of multivectors. See also [`lcontract`](@ref).

Equivalent to [`graded_prod(-, a, b)`](@ref).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ⨽ b` is the grade ``p - q`` part of `a*b`.
"""
rcontract(a, b) = graded_prod(-, a, b)

"$(@doc rcontract)"
⨽(a, b) = rcontract(a, b)



#= Exponentiation =#

# if a² is a scalar, then a²ⁿ = |a²|ⁿ and a²ⁿ⁺¹ = |a²|ⁿa
function power_with_scalar_square(a, a², p::Integer)
	# if p is even, p = 2n; if odd, p = 2n + 1
	aⁿ = a²^fld(p, 2)
	iseven(p) ? aⁿ*one(a) : aⁿ*a
end

function power_by_squaring(a::CompositeMultivector{Sig,S}, p::Integer) where {Sig,S}
	Π = one(Multivector{Sig,S})
	aⁿ = a
	while p > 0
		if isone(p & 1)
			Π *= aⁿ
		end
		aⁿ *= aⁿ
		p >>= 1
	end
	Π
end

Base.:^(a::BasisBlade, p::Integer) = power_with_scalar_square(a, scalar(a*a), p)
Base.literal_pow(::typeof(^), a::BasisBlade{Sig}, ::Val{2}) where {Sig} = BasisBlade{Sig,0}(0 => geometric_square_factor(Sig, bitsof(a))*a.coeff^2)

function Base.:^(a::CompositeMultivector{Sig,S}, p::Integer) where {Sig,S}
	# TODO: type stability?
	p < 0 && return inv(a)^abs(p)
	p == 0 && return one(a)
	p == 1 && return a
	a² = a*a
	p == 2 && return a²
	if isscalar(a²)
		power_with_scalar_square(a, scalar(a²), p)
	else
		power_by_squaring(a, p)
	end
end



#= Reversion and Dualities =#

"""
	graded_multiply(f, a::AbstractMultivector)

Multiply the grade `k` part of `a` by `f(k)`.
"""
graded_multiply(f, a::Scalar) = f(0)*a
graded_multiply(f, a::HomogeneousMultivector) = f(grade(a))*a
function graded_multiply(f, a::Multivector{Sig}) where Sig
	comps = copy(a.comps)
	dim = dimension(Sig)
	for k ∈ 0:dim
		comps[multivector_slice(Val(dim), Val(k))] *= f(k)
	end
	Multivector{Sig}(comps)
end

"""
	~a
	reversion(a::AbstractMultivector)

Reversion of a multivector.

Reversion is an anti-automorphism defined by reversing
the order of the geometric product: `~(a*b) == ~b * ~a`.
For a `k`-vector, `~k == reversion_sign(k) == (-1)^(k*(k - 1)/2)`.

| `mod(k, 4)` | [`reversion_sign(k)`](@ref) |
|:-----------:|:-------------------:|
| ``0`` | ``+1`` |
| ``1`` | ``+1`` |
| ``2`` | ``-1`` |
| ``3`` | ``-1`` |

"""
reversion(a) = graded_multiply(reversion_sign, a)

"$(@doc reversion)"
Base.:~(a::AbstractMultivector) = reversion(a)


"""
	involution(a)

Involute of a multivector.

Involution is an automorphism defined by reflecting through the origin:
for homogeneous multivectors, `involution(a) == (-1)^grade(a)*a`.
"""
involution(a) = graded_multiply(k -> (-1)^k, a)

"""
	a'ᶜ
	clifford_conj(a)

Clifford conjugate of a multivector.

Equivalent to `reversion(involution(a))`.
"""
clifford_conj(a) = graded_multiply(a) do k
	(-1)^k*reversion_sign(k)
end

"$(@doc clifford_conj)"
var"'ᶜ"(a) = clifford_conj(a)


"""
	flipdual(a)

A dual of a multivector, for when the overall sign isn’t important.

For a unit basis blade `a::BasisBlade`, the flipdual satisfies
`a*flipdual(a) == ±I`
where `±I` is the unit pseudoscalar or its negative.

Computing the `flipdual` is cheap, and is its own inverse: for a `BasisBlade`,
its bits are flipped, and for a `CompositeMultivector`, the components vector
is simply reversed.

The `flipdual` is _metric independent_, but depends on a choice of basis.
It differs from the Hodge and Poincaré duals by a per-grade scalar factor.
This makes it useful in projective geometry, where scalar factors are largely arbitrary.

See also [`hodgedual`](@ref) and [`poincaredual`](@ref).
"""
function flipdual end

"""
	poincaredual(a)

Poincaré dual of a multivector.

For a unit basis blade `a::BasisBlade`, the Poincaré dual satisfies
`a*poincaredual(a) == I`
where `I` is the unit pseudoscalar. I.e., `poincaredual(a)` is a “right complement” of `a`.

The Poincaré dual is _metric independent_, but depends on a choice of basis.
This makes is useful in degenerate algebras: non-zero multivectors have
non-zero Poincaré duals, even if their Hodge dual is zero.

See also [`hodgedual`](@ref) and [`flipdual`](@ref).

# Examples
```jldoctest
julia> @basis Cl(2,0,1)
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> hodgedual(v3), v3*v123 # Hodge dual is zero because v3*v3 == 0
(0, 0)

julia> poincaredual(v3)
BasisBlade{Cl(2,0,1), 2, Int64}:
 1 v12
```
"""
function poincaredual end


"""
	hodgedual(a) = ~a*I

Hodge dual of a multivector.

The Hodge dual is given by
```math
⋆a = ã I
```
where ``ã`` is the reversion of ``a`` and ``I`` is the unit pseudoscalar.
For ``k``-vectors ``a`` and ``b``, it is alternatively defined by
```math
a ∧ ⋆b = ⟨a, b⟩ I
```
where ``⟨a, b⟩ = a ⊙ b̃`` is the induced inner product on ``k``-vectors.

The Hodge dual depends on the _metric and orientation_ (choice of pseudoscalar).

See also [`poincaredual`](@ref).

# Examples
```jldoctest
julia> u = KVector{3,1}(1:3)
3-component KVector{3, 1, UnitRange{Int64}}:
 1 v1
 2 v2
 3 v3

julia> hodgedual(u)
3-component KVector{3, 2, Vector{Int64}}:
  3 v12
 -2 v13
  1 v23

julia> u ∧ hodgedual(u), u ⊙ ~u
(14v123, 14)

```
"""
function hodgedual end

for (name, signrule) in [
		:flipdual => (sig, bits) -> 1,
		:poincaredual => (sig, bits) -> sign_from_swaps(bits, bits_dual(dimension(sig), bits)),
		:hodgedual => (sig, bits) -> reversion_sign(count_ones(bits))geometric_prod_factor(sig, bits, first(bits_of_grade(dimension(sig)))),
	]
	@eval begin
		function $name(a::BasisBlade{Sig}) where {Sig}
			BasisBlade{signature(a)}(bits_dual(dimension(a), bitsof(a)) => $signrule(Sig, bitsof(a))*a.coeff)
		end

		function $name(a::KVector{Sig,K}) where {Sig,K}
			K′ = dimension(Sig) - K
			KVector{Sig,K′}(reverse($signrule.(Ref(Sig), bitsof(a)) .* a.comps))
		end

		function $name(a::Multivector{Sig}) where {Sig}
			Multivector{Sig}(reverse($signrule.(Ref(Sig), bitsof(a)) .* a.comps))
		end
	end
end




function outermorphism(mat::AbstractMatrix, a::AbstractMultivector{Sig}) where {Sig}
	resulttype = a isa HomogeneousMultivector ? KVector{Sig,grade(a)} : Multivector{Sig}
	a′ = zero(similar(resulttype, a))
	for (bits, coeff) ∈ nonzero_components(a)
		vs = KVector{Sig,1}.(eachcol(coeff*mat[:,bits_to_indices(bits)]))
		v = reduce(∧, vs; init = one(a′))
		add!(a′, v)
	end
	a′
end