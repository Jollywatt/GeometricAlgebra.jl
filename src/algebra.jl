#= Equality =#

Base.:(==)(a::BasisBlade{Sig}, b::BasisBlade{Sig}) where Sig = a.bits == b.bits ? a.coeff == b.coeff : iszero(a) && iszero(b)
Base.:(==)(a::Multivector{Sig,K}, b::Multivector{Sig,K}) where {Sig,K} = a.comps == b.comps
function Base.:(==)(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig}
	for k ∈ 0:dimension(a)
		if grade(a) ∋ k ∈ grade(b)
			a.comps[componentindices(a, k)] == b.comps[componentindices(b, k)] || return false
		elseif k ∈ grade(a)
			iszero(a.comps[componentindices(a, k)]) || return false
		elseif k ∈ grade(b)
			iszero(b.comps[componentindices(b, k)]) || return false
		end
	end
	true
end
Base.:(==)(a::BasisBlade{Sig}, b::Multivector{Sig}) where {Sig} = Multivector(a) == b
Base.:(==)(a::Multivector{Sig}, b::BasisBlade{Sig}) where {Sig} = a == Multivector(b)

Base.:(==)(a::AbstractMultivector, b::Scalar) = iszero(b) ? iszero(a) : isscalar(a) && scalar(a) == b
Base.:(==)(a::Scalar, b::AbstractMultivector) = iszero(a) ? iszero(b) : isscalar(b) && a == scalar(b)



#= Approximate Equality =#

isapproxzero(a; kwargs...) = isapprox(a, zero(a); kwargs...)
isapproxzero(a::BasisBlade; kwargs...) = isapproxzero(a.coeff; kwargs...)
isapproxzero(a::Multivector; kwargs...) = isapproxzero(a.comps; kwargs...)

Base.isapprox(a::BasisBlade{Sig}, b::BasisBlade{Sig}; kwargs...) where Sig = a.bits == b.bits ? isapprox(a.coeff, b.coeff; kwargs...) : isapproxzero(a) && isapproxzero(b)
function Base.isapprox(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}; kwargs...) where {Sig}
	k = promote_grades(a, b)
	isapprox(grade(a, k).comps, grade(b, k).comps; kwargs...)
end

Base.:isapprox(a::AbstractMultivector, b::Scalar; kwargs...) = isapprox(a, zero(a) + b; kwargs...)
Base.:isapprox(a::Scalar, b::AbstractMultivector; kwargs...) = isapprox(zero(b) + a, b; kwargs...)



#= Scalar Multiplication =#

scalar_multiply(a::BasisBlade{Sig,K}, b) where {Sig,K} = BasisBlade{Sig,K}(a.coeff*b, a.bits)
scalar_multiply(a, b::BasisBlade{Sig,K}) where {Sig,K} = BasisBlade{Sig,K}(a*b.coeff, b.bits)

scalar_multiply(a::Multivector, b) = constructor(a)(a.comps.*b)
scalar_multiply(a, b::Multivector) = constructor(b)(a.*b.comps)

Base.:*(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
Base.:*(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)
Base.:-(a::AbstractMultivector) = -numberone(eltype(a))*a

promote_to(T, x) = convert(promote_type(T, typeof(x)), x)
Base.:/(a::AbstractMultivector, b::Scalar) = a*inv(promote_to(eltype(a), b))
Base.:\(a::Scalar, b::AbstractMultivector) = inv(promote_to(eltype(b), a))*b

Base.://(a::AbstractMultivector, b::Scalar) = a*(one(b)//b)
Base.://(a::Scalar, b::AbstractMultivector) = inv(b//a)
Base.://(a::AbstractMultivector, b::AbstractMultivector) = a*(numberone(eltype(b))//b)



#= Addition =#

"""
	add!(a::Multivector, b::Blade)
	add!(a::Multivector, bits, coeff)

Add the blade coefficient to the corresponding component of a multivector,
if the multivector has such a component.

!!! warning
	If the multivector cannot represent components of the required grade, it is returned unmodified.

This mutates and returns `a` if it is a mutable type, otherwise it returns a new multivector of identical type.
(Thus, the blade coefficient must be convertible to the multivector’s eltype.)
"""
function add!(a::Multivector, coeff, bits::Unsigned)
	i = componentindex(a, bits)
	isnothing(i) && return a
	if issetindexable(a.comps)
		a.comps[i] += coeff
		a
	else
		comps = setindex(a.comps, a.comps[i] + coeff, i)
		constructor(a)(comps)
	end
end

add!(a::Multivector, b::BasisBlade) = add!(a, b.coeff, b.bits)

function add!(a::Multivector{Sig,K}, b::Multivector{Sig,K}) where {Sig, K}
	if issetindexable(a.comps)
		a.comps .+= b.comps
		a
	else
		typeof(a)(a.comps + b.comps)
	end
end

function add!(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig}
	if issetindexable(a.comps)
		for k ∈ 0:dimension(Sig)
			grade(a) ∋ k ∈ grade(b) || continue
			a.comps[componentindices(a, k)] .+= b.comps[componentindices(b, k)]
		end
	else
		for (coeff, bits) in nonzero_components(b)
			a = add!(a, coeff, bits)
		end
	end
	a
end

add!(Σ, a, b, c...) = add!(add!(Σ, a), b, c...)

resulting_grades(::typeof(+), dim, pq...) = promote_grades(dim, pq...)

Base.:+(a::Multivector{Sig,K}, b::Multivector{Sig,K}) where {Sig,K} = Multivector{Sig,K}(a.comps + b.comps)
function Base.:+(a::AbstractMultivector{Sig}, bc::AbstractMultivector{Sig}...) where {Sig}
	Σ = zero(resulting_multivector_type(+, a, bc...))
	add!(Σ, a, bc...)
end
Base.:-(a::Multivector{Sig,K}, b::Multivector{Sig,K}) where {Sig,K} = Multivector{Sig,K}(a.comps - b.comps)
Base.:-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)



#= Scalar Addition =#

add_scalar(a::AbstractMultivector{Sig}, b::Scalar) where {Sig} = a + BasisBlade{Sig}(b)

Base.:+(a::AbstractMultivector, b::Scalar) = add_scalar(a, b)
Base.:+(a::Scalar, b::AbstractMultivector) = add_scalar(b, a)

Base.:-(a::AbstractMultivector, b::Scalar) = add_scalar(a, -b)
Base.:-(a::Scalar, b::AbstractMultivector) = add_scalar(-b, a)




#= Geometric Multiplication =#

"""
	a * b
	geometric_prod(a, b)

Geometric product of multivectors.
"""
geometric_prod(a::Scalar, b::Scalar) = a*b
geometric_prod(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
geometric_prod(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)

function geometric_prod(a::BasisBlade{Sig}, b::BasisBlade{Sig}) where {Sig}
	factor = geometric_prod_factor(Sig, a.bits, b.bits)
	BasisBlade{Sig}(factor*(a.coeff*b.coeff), a.bits ⊻ b.bits)
end

function resulting_grades(::typeof(geometric_prod), dim, p::Integer, q::Integer)
	dim == max(p, q) && return dim - min(p, q)
	abs(p - q):2:min(p + q, dim)
end

@symbolic_optim function geometric_prod(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	c = zero(resulting_multivector_type(geometric_prod, a, b))
	for (acoeff, abits) ∈ nonzero_components(a), (bcoeff, bbits) ∈ nonzero_components(b)
		# println(abits, ":", bbits)
		factor = geometric_prod_factor(Sig, abits, bbits)
		c = add!(c, factor*(acoeff*bcoeff), abits ⊻ bbits)
	end
	c
end

Base.:*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)



#= Scalar Product =#

"""
	a ⊙ b
	scalar_prod(a, b) -> Number

Scalar part of the geometric product `a*b`.
"""
scalar_prod(a::Scalar, b::Scalar) = a*b
scalar_prod(a::AbstractMultivector, b::Scalar) = scalar(a)*b
scalar_prod(a::Scalar, b::AbstractMultivector) = a*scalar(b)

function scalar_prod(a::BasisBlade{Sig,K}, b::BasisBlade{Sig,K}) where {Sig,K}
	if a.bits == b.bits
		geometric_prod_factor(Sig, a.bits, b.bits)*(a.coeff*b.coeff)
	else
		numberzero(promote_type(eltype(a), eltype(b)))
	end
end
scalar_prod(a::BasisBlade{Sig}, b::BasisBlade{Sig}) where {Sig} = numberzero(promote_type(eltype(a), eltype(b)))

function scalar_prod(a::Multivector{Sig,K}, b::Multivector{Sig,K}) where {Sig,K}
	s = numberzero(promote_type(eltype(a), eltype(b)))
	for (acoeff, bcoeff, bits) in zip(a.comps, b.comps, componentbits(Multivector{Sig,K}))
		s += geometric_square_factor(Sig, bits)*(acoeff*bcoeff)
	end
	s
end

@symbolic_optim function scalar_prod(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	s = numberzero(promote_type(eltype(a), eltype(b)))
	for k in 0:dimension(Sig)
		grade(a) ∋ k ∈ grade(b) || continue
		s += scalar_prod(grade(a, k), grade(b, k))
	end
	s
end

@doc (@doc scalar_prod)
a ⊙ b = scalar_prod(a, b)



#= Graded Products =#

"""
	graded_prod(grade_selector::Function, a, b)

A grade-filtered product of multivectors.

Terms of the geometric product `a*b` are filtered according to `grade_selector`.
For instance, if `grade(a) == p` and `grade(b) == q`, then `graded_prod(f, a, b)` is the grade `f(p, q)` part of `a*b`.
The extends linearly to general multivectors ``A`` and ``B`` as
```math
	(A, B) ↦ \\sum_{p,q} ⟨⟨A⟩_p ⟨B⟩_q⟩_{f(p, q)}
```
where ``⟨⋅⟩_k`` denotes the grade ``k`` part.

The [wedge ``∧``](@ref GeometricAlgebra.wedge), [inner ``⋅``](@ref GeometricAlgebra.inner) and
[contraction](@ref GeometricAlgebra.lcontract) products are special cases of this product.
For example, the wedge product is defined as:
```julia
wedge(a, b) = graded_prod(+, a, b)
```
"""
graded_prod(grade_selector, a::Scalar, b::Scalar) = iszero(grade_selector(0, 0)) ? a*b : zero(promote_type(a, b))
graded_prod(grade_selector, a::AbstractMultivector, b::Scalar) = grade(a) == (grade_selector(grade(a), 0)) ? scalar_multiply(a, b) : zero(a)
graded_prod(grade_selector, a::Scalar, b::AbstractMultivector) = grade(b) == (grade_selector(0, grade(b))) ? scalar_multiply(a, b) : zero(a)

function graded_prod(grade_selector::Function, a::BasisBlade{Sig}, b::BasisBlade{Sig}) where {Sig}
	K = grade_selector(grade(a), grade(b))
	T = promote_type(eltype(a), eltype(b))
	0 <= K <= dimension(Sig) || return BasisBlade{Sig}(numberzero(T))

	if count_ones(a.bits ⊻ b.bits) == K
		BasisBlade{Sig,K,T}(geometric_prod_factor(Sig, a.bits, b.bits)*(a.coeff*b.coeff), a.bits ⊻ b.bits)
	else
		bits = UInt(1) << K - 1 # the coeff is zero, so these are arbitrary — but for consistency they should match the grade
		BasisBlade{Sig,K,T}(numberzero(T), bits)
	end
end

resulting_grades(::Tuple{typeof(graded_prod),GradeSelector}, dim, p::Integer, q::Integer) where {GradeSelector} = GradeSelector.instance(p, q)

@symbolic_optim function graded_prod(grade_selector::Function, a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	c = zero(resulting_multivector_type((graded_prod, grade_selector), a, b))
	for (acoeff, abits) ∈ nonzero_components(a), (bcoeff, bbits) ∈ nonzero_components(b)
		bits = abits ⊻ bbits
		if count_ones(bits) == grade_selector(count_ones(abits), count_ones(bbits))
			factor = geometric_prod_factor(Sig, abits, bbits)
			c = add!(c, factor*(acoeff*bcoeff), bits)
		end
	end
	c
end

@symbolic_optim function graded_prod(::Val{K}, a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig,K}
	c = zero(Multivector{Sig,K}, promote_type(eltype(a), eltype(b)))
	for (acoeff, abits) ∈ nonzero_components(a), (bcoeff, bbits) ∈ nonzero_components(b)
		bits = abits ⊻ bbits
		if count_ones(bits) in K
			factor = geometric_prod_factor(Sig, abits, bbits)
			c = add!(c, factor*(acoeff*bcoeff), bits)
		end
	end
	c
end


#= Derived Products =#

"""
	a ∧ b
	wedge(a, b, ...)

Wedge product of multivectors (a.k.a. the _outer_, _exterior_, _progressive_ or _alternating_ product, or _join_).

This is a grade-raising operation, equivalent to [`graded_prod(+, a, b)`](@ref GeometricAlgebra.graded_prod).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ∧ b` is defined as
the grade ``p + q`` part of `a*b`.
"""
wedge(a, b) = graded_prod(+, a, b)
@doc (@doc wedge)
a ∧ b = wedge(a, b)

# since wedge is associative, this can be handy
wedge(a, b, c, d...) = wedge(wedge(a, b), c, d...)
wedge(a) = a

"""
	a ⋅ b
	inner(a, b)

Inner product of multivectors.

This is a grade lowering operation, equivalent to [`graded_prod(abs∘-, a, b)`](@ref GeometricAlgebra.graded_prod).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ⋅ b` is defined as
the grade ``|p - q|`` part of `a*b`.

Note that for scalars `a` and `b`, the inner product reduces to scalar multiplication,
in contrast to some authors (see [Dorst2002](@cite) for discussion).

See also [`lcontract`](@ref) and [`rcontract`](@ref).
"""
inner(a, b) = graded_prod(abs∘-, a, b)
@doc (@doc inner)
⋅(a, b) = inner(a, b)


"""
	a ⨼ b
	lcontract(a, b)

Left contraction of multivectors.

Equivalent to [`graded_prod((p, q) -> q - p, a, b)`](@ref GeometricAlgebra.graded_prod).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ⨼ b` is defined as
the grade ``q - p`` part of `a*b`.

See also [`rcontract`](@ref) and [`inner`](@ref).
"""
lcontract(a, b) = graded_prod((-)∘-, a, b)
@doc (@doc lcontract)
⨼(a, b) = lcontract(a, b)


"""
	a ⨽ b
	rcontract(a, b)

Right contraction of multivectors.

Equivalent to [`graded_prod(-, a, b)`](@ref GeometricAlgebra.graded_prod).
If `a` and `b` are of grades ``p`` and ``q`` respectively, then `a ⨽ b` is defined as
the grade ``p - q`` part of `a*b`.

See also [`lcontract`](@ref) and [`inner`](@ref).
"""
rcontract(a, b) = graded_prod(-, a, b)
@doc (@doc rcontract)
⨽(a, b) = rcontract(a, b)



#= Exponentiation =#

# if a² is a scalar, then a²ⁿ = |a²|ⁿ and a²ⁿ⁺¹ = |a²|ⁿa
function power_with_scalar_square(a, a², p::Integer)
	# if p is even, p = 2n; if odd, p = 2n + 1
	aⁿ = a²^fld(p, 2)
	iseven(p) ? aⁿ*one(a) : aⁿ*a
end

square(a::BasisBlade{Sig}) where Sig = BasisBlade{Sig}(geometric_square_factor(Sig, a.bits)*a.coeff^2)
function square(a::Multivector)
	dim = dimension(a)
	# (pseudo)(scalars|vectors) always square to scalars
	k = grade(a)
	if k == 0 || k == 1 || k == dim || k == dim - 1
		inner(a, a) # we want this branch to be taken at compile time
	else 
		a*a
	end
	# grade(a) ∈ (0, 1, dim - 1, dim) ? inner(a, a) : a*a
end

Base.literal_pow(::typeof(^), a::AbstractMultivector, ::Val{2}) = square(a)

Base.:^(a::BasisBlade, p::Integer) = power_with_scalar_square(a, scalar(a*a), p)
function Base.:^(a::Multivector{Sig,S}, p::Integer) where {Sig,S}
	p < 0 && return inv(a)^abs(p)
	p == 0 && return one(a)
	p == 1 && return a
	a² = square(a)
	p == 2 && return a²
	if isscalar(a²)
		power_with_scalar_square(a, scalar(a²), p)
	else
		Base.power_by_squaring(a, p)
	end
end



#= Reversion and Dualities =#

"""
	graded_multiply(f, a::AbstractMultivector)

Multiply the grade `k` part of `a` by `f(k)`.
"""
graded_multiply(f, a::Scalar) = f(0)a
function graded_multiply(f, a::AbstractMultivector)
	if ishomogeneous(a)
		f(grade(a))a
	else
		constructor(a)(f.(count_ones.(componentbits(a))) .* a.comps)
	end
end

"""
	~a
	reversion(a::AbstractMultivector)

Reversion of a multivector.

Reversion is an anti-automorphism defined by reversing
the order of the geometric product: `~(a*b) == ~b * ~a`.
For a `k`-vector `a`, the reversion is `reversion_sign(k)*a`
where the sign is given by ``(-1)^{k(k - 1)/2}``.

See also [`involution`](@ref) and [`clifford_conj`](@ref).
"""
reversion(a) = graded_multiply(reversion_sign, a)
@doc (@doc reversion)
Base.:~(a::AbstractMultivector) = reversion(a)


"""
	involution(a)

Involute of a multivector.

Involution is an automorphism defined by reflecting through the origin:
for homogeneous multivectors, `involution(a) == (-1)^grade(a)*a`.

See also [`reversion`](@ref) and [`clifford_conj`](@ref).
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
@doc (@doc clifford_conj)
var"'ᶜ"(a) = clifford_conj(a)


"""
	flipdual(a)

A dual of a multivector, for when the overall sign isn’t important.

For a unit `a::BasisBlade`, the flipdual satisfies
`a*flipdual(a) == ±I`
where `±I` is the unit pseudoscalar or its negative.

The `flipdual` is cheap to compute and is its own inverse.
It simply flips the bits of a `BasisBlade`, or reverses the components
vector of a `Multivector`.

The `flipdual` is _metric independent_ (but depends on a choice of _orientation_, or the order of basis vectors).

See also [`hodgedual`](@ref).
"""
function flipdual end


"""
	ldual(a)
	rdual(a)

Left and right multivector duals (a.k.a., _complements_).
The right dual is also called the _Poincaré dual_.

For a unit basis blade `a`, the duals satisfy `a*rdual(a) == I == ldual(a)*a` where `I` is the unit pseudoscalar.
If `dimension(a)` is odd, `rdual` and `ldual` are identical and self-inverse; in general, they are inverses of each other.

The left and right duals are _metric independent_ (but depend on a choice of _orientation_, or the order of basis vectors).
This makes them useful in degenerate algebras where `I^2 == 0`, since a non-zero multivector always
has a non-zero dual, even if its Hodge dual is zero.

See also [`hodgedual`](@ref).
"""
function ldual end

@doc (@doc ldual)
function rdual end



"""
	hodgedual(a) = ~a*I

Hodge dual of a multivector.

The Hodge dual is defined by
```math
H(a) = ã I
```
where ``ã`` is the reversion of ``a`` and ``I`` is the unit pseudoscalar.
For ``k``-vectors ``a`` and ``b``, it is alternatively defined by
```math
a ∧ H(b) = ⟨a, b⟩ I
```
where ``⟨a, b⟩ = a ⊙ b̃`` is the induced inner product on ``k``-vectors.

The Hodge dual is _metric dependent_, since it involves multiplication by `I`.

See also [`invhodgedual`](@ref) and [`ldual`](@ref), [`rdual`](@ref).

# Examples
```jldoctest
julia> u = Multivector{3,1}(1:3)
3-component Multivector{3, 1, UnitRange{Int64}}:
 1 v1
 2 v2
 3 v3

julia> hodgedual(u)
3-component Multivector{3, 2, SVector{3, Int64}}:
  3 v12
 -2 v13
  1 v23

```
"""
function hodgedual end

"""
	invhodgedual(a)

Inverse of the multivector Hodge dual.

In degenerate algebras (for which ``I^2 = 0``), the Hodge dual is not invertible.
However, if `a` is a basis blade with a non-zero Hodge dual, then `invhodgedual(hodgedual(a)) == a` holds.

See also [`hodgedual`](@ref).
"""
function invhodgedual end

for (name, signrule) in [
		:flipdual => (sig, bits) -> 1,
		:rdual => (sig, bits) -> sign_from_swaps(bits, bits_dual(dimension(sig), bits)),
		:ldual => (sig, bits) -> sign_from_swaps(bits_dual(dimension(sig), bits), bits),
		:hodgedual => function(sig, bits)
			I = first(bits_of_grade(dimension(sig)))
			reversion_sign(count_ones(bits))geometric_prod_factor(sig, bits, I)
		end,
		:invhodgedual => function(sig, bits)
			I = first(bits_of_grade(dimension(sig)))
			reversion_sign(count_ones(I ⊻ bits))geometric_prod_factor(sig, I ⊻ bits, I)
		end,
	]
	@eval begin
		function $name(a::BasisBlade{Sig}) where {Sig}
			bits = bits_dual(dimension(a), a.bits)
			BasisBlade{signature(a)}($signrule(Sig, a.bits)*a.coeff, bits)
		end

		@symbolic_optim function $name(a::Multivector{Sig,K}) where {Sig,K}
			K′ = promote_grades(dimension(Sig), dimension(Sig) .- K)
			Multivector{Sig,K′}(reverse($signrule.(Ref(Sig), componentbits(a)) .* a.comps))
		end
	end
end


"""
	a ∨ b
	antiwedge(a, b)

Anti-wedge product of multivectors (a.k.a. the _regressive_ product or _meet_).

The anti-wedge product satisfies
```math
D(a ∨ b) = (D a) ∧ (D b)
```
where ``D`` is a duality operation such as `ldual`, `rdual` or, if ``I^2 ≠ 0``, `hodgedual`.

The anti-wedge product is _metric independent_ like the wedge product,
but does depend on the choice of _orientation_ (the ordering of basis vectors).
"""
antiwedge(a, b) = rdual(ldual(a)∧ldual(b))

@doc (@doc antiwedge)
∨(a, b) = antiwedge(a, b)


"""
	sandwich_prod(R, a)

Sandwich product `R*a*~R` of multivector `a` by a rotor `R`.
"""
function sandwich_prod end

sandwich_prod(R, a) = grade(R*a*reversion(R), grade(a))
@symbolic_optim sandwich_prod(R::AbstractMultivector, a::AbstractMultivector) = grade(R*a*reversion(R), grade(a))


"""
	outermorphism(mat, a)

Outermorphism of the multivector `a` specified by the matrix `mat`.

If ``f`` is a linear map, then the outermorphism ``f̲`` is a linear map  satisfying
``f̲(𝒖) = f(𝒖)`` on vectors ``𝒖`` and ``f̲(a ∧ b) = f̲(a) ∧ f̲(b)`` on
general multivectors.
"""
function outermorphism(mat::AbstractMatrix, a::AbstractMultivector{Sig}; sig=Sig) where {Sig}
	a′ = zero(Multivector{sig,grade(a)}, promote_type(eltype(mat), eltype(a)))
	for (coeff, bits) ∈ nonzero_components(a)
		vs = Multivector{sig,1}.(eachcol(mat[:,bits_to_indices(bits)]))
		v = reduce(∧, vs; init = one(a′))
		a′ = add!(a′, coeff*v)
	end
	a′
end
outermorphism(mat, a::Scalar) = a

"""
	embed(sig, a::Multivector)

Embed a multivector into the algebra of metric signature `sig`, adding or discarding dimensions as necessary.

Basis vectors in the original and new spaces are associated by the order in which they appear.
Components are dropped if `dimension(sig) < dimension(a)`, and extra zero components are added
if `dimension(sig) > dimension(a)`.

# Examples
```jldoctest
julia> embed(Cl(2,2), Multivector{3,2}([1, 2, 3]))
6-component Multivector{Cl(2,2), 2, SVector{6, Float64}}:
 1.0 v12
 2.0 v13
 3.0 v23
 0.0 v14
 0.0 v24
 0.0 v34

julia> embed(2, Multivector{3,1}([1,2,3]))
2-component Multivector{2, 1, SVector{2, Int64}}:
 1 v1
 2 v2
```
"""
function embed(sig, a::Multivector)
	T = Multivector{sig,grade(a)}
	b = zero(T)
	for k in grade(a)
		n_orig = ncomponents(signature(a), k)
		n_new = ncomponents(sig, k)
		if n_new > n_orig
			b += Multivector{sig,k}([a[k].comps; zeros(n_new - n_orig)])
		else
			b += Multivector{sig,k}(a[k].comps[1:n_new])
		end
	end
	b
end
embed(sig, a::BasisBlade) = BasisBlade{sig}(iszero(a.bits >> dimension(sig)) ? a.coeff : zero(a.coeff), a.bits)

unit_pseudoscalar(::Val{Sig}) where {Sig} = let dim = dimension(Sig)
	BasisBlade{Sig,dim}(1, bits_dual(dim, UInt(0)))
end
unit_pseudoscalar(a::AbstractMultivector{Sig}) where {Sig} = unit_pseudoscalar(Val(Sig))

pseudoscalar_square(sig::Val) = scalar(unit_pseudoscalar(sig)^2)
pseudoscalar_square(a::AbstractMultivector{Sig}) where {Sig} = pseudoscalar_square(Val(Sig))

