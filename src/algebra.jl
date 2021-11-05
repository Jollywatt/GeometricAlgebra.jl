#= EQUALITY

Only equality between multivectors of the *almost* the same type is defined,
leveraging the promotion interface to enable comparison between differing types.

Multivector promotion (intentionally) does *not* unify the grade and bits parameters,
so the following methods must cover types which are identical except for those parameters.

Note that two blades of differing basis or two multivectors of differing grade are
treated as equal if they are both zero.
=#

==(a::Blade{sig}, b::Blade{sig}) where sig = bitsof(a) == bitsof(b) ? a.coeff == b.coeff : iszero(a) && iszero(b)
==(a::(Multivector{sig,k,S} where k), b::(Multivector{sig,k,S} where k)) where {sig,S} = grade(a) == grade(b) ? a.components == b.components : iszero(a) && iszero(b)
==(a::MixedMultivector{sig,S}, b::MixedMultivector{sig,S}) where {sig,S} = a.components == b.components

==(a::AbstractMultivector, b::Scalar) = iszero(b) && iszero(a) || isscalar(a) && scalar(a) == b
==(a::Scalar, b::AbstractMultivector) = iszero(a) && iszero(b) || isscalar(b) && a == scalar(b)

==(a::AbstractMultivector{sig}...) where sig = ==(promote(a...)...)
==(a::AbstractMultivector...) = false # multivectors with non-identical signatures are teated as non-equal


function Base.isapprox(a::Blade{sig}, b::Blade{sig}; kwargs...) where sig
	if bitsof(a) == bitsof(b)
		isapprox(a.coeff, b.coeff; kwargs...)
	else
		isapprox(a.coeff, zero(a.coeff); kwargs...) && isapprox(b.coeff, zero(a.coeff); kwargs...)
	end
end
Base.isapprox(a::T, b::T; kwargs...) where {T<:MixedMultivector} = isapprox(a.components, b.components; kwargs...)
function Base.isapprox(a::Multivector{sig}, b::Multivector{sig}; kwargs...) where sig
	if grade(a) == grade(b)
		isapprox(a.components, b.components; kwargs...)
	else
		# multivectors of different grade are approximately equal is they are both approximately zero
		isapprox(a, zero(a); kwargs...) && isapprox(b, zero(b); kwargs...)
	end
end
Base.isapprox(a::AbstractMultivector, b::AbstractMultivector; kwargs...) = isapprox(promote(a, b)...; kwargs...)

Base.isapprox(a::AbstractMultivector, b::Scalar; kwargs...) = isapprox(a, one(a)*b; kwargs...)
Base.isapprox(a::Scalar, b::AbstractMultivector; kwargs...) = isapprox(a*one(b), b; kwargs...)






#= SCALAR MULTIPLICATION =#

scalar_multiply(a::Blade{sig}, b) where {sig} = Blade{sig}(a.coeff*b, bitsof(a))
scalar_multiply(a, b::Blade{sig}) where {sig} = Blade{sig}(a*b.coeff, bitsof(b))

constructor(::Multivector{sig,k}) where {sig,k} = Multivector{sig,k}
constructor(::MixedMultivector{sig}) where sig = MixedMultivector{sig}
scalar_multiply(a::CompositeMultivector{<:AbstractVector}, b) = constructor(a)(a.components*b)
scalar_multiply(a, b::CompositeMultivector{<:AbstractVector}) = constructor(b)(a*b.components)

*(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
*(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)

-(a::AbstractMultivector) = -one(eltype(a))*a

promote_to(x::T, ::Type{S}) where {T,S} = convert(promote_type(T, S), x)
/(a::AbstractMultivector, b::Scalar) = a*inv(promote_to(b, eltype(a)))
\(a::Scalar, b::AbstractMultivector) = inv(promote_to(a, eltype(b)))*b



#= ADDITION =#


scalar_add(a::Blade{sig,0}, b) where {sig} = Blade{sig}(a.coeff + b, bitsof(a))
scalar_add(a::Multivector{sig,0}, b) where {sig} = Multivector{sig,0}(a.components .+ b)
scalar_add(a::HomogeneousMultivector, b) = scalar_add(MixedMultivector(a), b)

scalar_add!(a::MixedMultivector, b) = (a.components[begin] += b; a)
function scalar_add(a::MixedMultivector{sig,<:AbstractVector}, b) where {sig}
	# this allows eltype promotion (unlike, e.g., copying and setting first element)
	comps = [a.components[begin] + b; a.components[begin + 1:end]]
	best_type(a, set_eltype=eltype(comps))(comps)
end




# warning: `for u ∈ blades(b) ... end` is not fast

# mutating addition which does not require same element/storage types
function add!(a::CompositeMultivector, b::Blade)
	a.components[bits_to_index(a, bitsof(b))] += b.coeff
	a
end
function add!(a::MixedMultivector, b::CompositeMultivector)
	a.components += full_components_vector(b)
	a
end

# mutating add! falls back to normal add for immutable storagetypes
add!(a::CompositeMultivector{<:StaticVector}, b::Blade) = add(a, b)
add!(a::CompositeMultivector{<:StaticVector}, b::Multivector) = add(a, b)
add!(a::CompositeMultivector{<:StaticVector}, b::MixedMultivector) = add(a, b)

function add(As::Multivector{sig,k,<:AbstractVector}...) where {sig,k}
	Multivector{sig,k}(.+((a.components for a ∈ As)...))
end
function add(As::MixedMultivector{sig,<:AbstractVector}...) where sig
	MixedMultivector{sig}(.+((a.components for a ∈ As)...))
end

# dispatches on blades/multivectors with uniform grade
function add(As::HomogeneousMultivector{sig,k}...) where {sig,k}
	add(convert.(best_type(Multivector, As...), As)...)
end
function add(As::AbstractMultivector...)
	add(convert.(best_type(MixedMultivector, As...), As)...)
end



+(a::AbstractMultivector, b::Scalar) = scalar_add(a, b)
+(a::Scalar, b::AbstractMultivector) = scalar_add(b, a)
+(a::AbstractMultivector) = a
+(as::AbstractMultivector...) = add(as...)

-(a::AbstractMultivector, b::Scalar) = a + (-b)
-(a::Scalar, b::AbstractMultivector) = a + (-b)
-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)






#= THE GEOMETRIC PRODUCT =#

"""
	a*b
	geometric_prod(a, b)

Geometric product of multivectors.

If both arguments are blades, then the result is a `Blade`.
Otherwise, the geometric product is general `MixedMultivector`.
"""
geometric_prod

function geometric_prod(a::Blade, b::Blade)
	sig = shared_sig(a, b)
	factor, bits = geometric_prod_bits(sig, bitsof(a), bitsof(b))
	Blade{sig}(factor*(a.coeff*b.coeff), bits)
end

function geometric_prod(a::AbstractMultivector, b::AbstractMultivector)
	Π = zero(best_type(MixedMultivector, a, b))
	for (abits, acomp) ∈ _blades(a), (bbits, bcomp) ∈ _blades(b)
		factor, bits = geometric_prod_bits(signature(Π), abits, bbits)
		Π.components[bits_to_mmv_index(bits, dimension(Π))] += factor*acomp*bcomp
	end
	Π
end


"""
	*(::AbstractMultivector, ::AbstractMultivector)

Geometric product. See [`$(repr(geometric_prod))`](@ref).
"""
*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)



#= GRADED PRODUCTS =#

"""
	homogeneous_prod(a, b, Val(k))

The grade `k` part of the geometric product of `a` and `b`, i.e., ``⟨ab⟩_k``.
Equivalent to, but more efficient than, `grade(geometric_prod(a, b), k)`.
"""
homogeneous_prod

homogeneous_prod(a::Scalar, b::Scalar, ::Val{k}) where k = iszero(k) ? a*b : zero(a*b)
function homogeneous_prod(a::Blade, b::Blade, ::Val{k}) where k
	bits = bitsof(a) ⊻ bitsof(b)
	if grade(bits) == k
		geometric_prod(a, b)
	else
		zero(best_type(Blade, a, b; grade=Val(k)))
	end
end
function homogeneous_prod(a::AbstractMultivector, b::AbstractMultivector, ::Val{k}) where k
	Π = zero(best_type(Multivector, a, b; grade=Val(k)))

	if (a, b) isa NTuple{2,HomogeneousMultivector}
		k ∈ abs(grade(a) - grade(b)):2:(grade(a) + grade(b)) || return Π
	else
		0 <= k <= dimension(a) || return Π
	end

	for u ∈ blades(a), v ∈ blades(b)
		if grade(bitsof(u) ⊻ bitsof(v)) == k
			Π = add!(Π, homogeneous_prod(u, v, k))
		end
	end
	Π
end
homogeneous_prod(a, b, k::Integer) = homogeneous_prod(a, b, Val(k))

"""
	graded_prod(a, b, grade_selector::Function)

Homogenised geometric product with resulting grade given as a function of the grades of `a` and `b`.

If ``a_p`` and ``b_q`` are multivectors of grade ``p`` and ``q`` respectively, then this is
equal to ``⟨a_p b_q⟩_{k(p,q)}`` where ``k(p,q)`` is `grade_selector(p, q)`.
If `a` or `b` are multi-graded, then this extends linearly to ``∑_{p,q} ⟨⟨a⟩_p ⟨b⟩_q⟩_{k(p,q)}``.
"""
graded_prod

function graded_prod(a::AbstractMultivector, b::AbstractMultivector, grade_selector::Function)
	Π = zero(best_type(MixedMultivector, a, b))
	for u ∈ blades(a), v ∈ blades(b)
		k = grade_selector(grade(u), grade(v))
		0 <= k <= dimension(Π) || continue
		Π = add!(Π, homogeneous_prod(u, v, k))
	end
	Π
end
function graded_prod(a::HomogeneousMultivector, b::HomogeneousMultivector, grade_selector::Function)
	homogeneous_prod(a, b, grade_selector(grade(a), grade(b)))
end
graded_prod(a::AbstractMultivector, b::Scalar, grade_selector) = graded_prod(a, one(a)b, grade_selector)
graded_prod(a::Scalar, b::AbstractMultivector, grade_selector) = graded_prod(one(b)a, b, grade_selector)
graded_prod(a::Scalar, b::Scalar, grade_selector) = homogeneous_prod(a, b, grade_selector(0, 0))


"""
	∧(a, b)

Wedge product of multivectors.

For homogeneous multivectors, `a∧b` is of grade `grade(a) + grade(b)`.
For mixed multivectors, the definition extends by linearity to ``a∧b = ∑_{pq} ⟨⟨a⟩_p ⟨b⟩_q⟩_{p + q}``.
"""
∧(a, b) = graded_prod(a, b, +)

"""
	⋅(a, b)

Generalized inner product (or "fat" dot product) of multivectors.

For homogeneous multivectors, `a⋅b` is of grade `abs(grade(a) - grade(b))`.
For mixed multivectors, the definition extends by linearity to ``a⋅b = ∑_{pq} ⟨⟨a⟩_p ⟨b⟩_q⟩_{|p - q|}``.
"""
⋅(a, b) = graded_prod(a, b, abs∘-)

"""
	⨼(a, b)

Left contraction of multivectors.

For homogeneous multivectors, `a⨼b` is of grade `grade(b) - grade(a)`.
For mixed multivectors, the definition extends by linearity to ``a⨼b = ∑_{pq} ⟨⟨a⟩_p ⟨b⟩_q⟩_{q - p}``.
"""
⨼(a, b) = graded_prod(a, b, (p, q) -> q - p)

"""
	⨽(a, b)

Right contraction of multivectors.

For homogeneous multivectors, `a⨽b` is of grade `grade(a) - grade(b)`.
For mixed multivectors, the definition extends by linearity to ``a⨽b = ∑_{pq} ⟨⟨a⟩_p ⟨b⟩_q⟩_{p - q}``.
"""
⨽(a, b) = graded_prod(a, b, (p, q) -> p - q)

"""
	∗(a, b)

Scalar product of multivectors, equal to `scalar(a*b)`.
"""
∗(a, b) = graded_prod(a, b, (p, q) -> 0)




#= MULTIPLICATIVE INVERSES =#

Base.inv(a::Blade) = a/scalar(a*a) # blades are guaranteed to have scalar squares
function Base.inv(a::CompositeMultivector)
	a² = a^2
	if isscalar(a²)
		a/scalar(a²)
	else
		inv_matrixmethod(a)
	end
end

function linear_map(a::AbstractMultivector)
	n = dimension(a)
	hcat([full_components_vector(a*blade_like(a, 1, mmv_index_to_bits(i, n)))
		for i ∈ 1:2^dimension(a)]...)
end

function inv_matrixmethod(a::AbstractMultivector)
	n = dimension(a)
	A = hcat([full_components_vector(a*blade_like(a, 1, mmv_index_to_bits(i, n)))
		for i ∈ 1:2^dimension(a)]...)
	id = full_components_vector(one(a))
	inv_components = A\id
	# result should be converted into the same storage type as original
	a⁻¹ = best_type(MixedMultivector, a; set_eltype=eltype(inv_components))(inv_components)
end


/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b

/(a::Scalar, b::AbstractMultivector) = a*inv(b)
\(a::AbstractMultivector, b::Scalar) = inv(a)*b



#= EXPONENTIATION =#

# if a² is a scalar, then a²ⁿ = |a²|ⁿ, a²ⁿ⁺¹ = |a²|ⁿa
function power_with_scalar_square(a, a², p::Integer)
	# if p is even, p = 2n; if odd, p = 2n + 1
	aⁿ = a²^fld(p, 2)
	iseven(p) ? aⁿ*one(a) : aⁿ*a
end

function power_by_squaring(a::AbstractMultivector, p::Integer)
	p >= 0 || return power_by_squaring(inv(a), abs(p))
	Π = one(best_type(MixedMultivector, a))
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

^(a::Blade, p::Integer) = power_with_scalar_square(a, scalar(a*a), p)

function ^(a::CompositeMultivector, p::Integer)
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



#= GRADE PROJECTION =#

grade(a::Scalar) = 0
grade(a, k::Integer) where sig = grade(a) == k ? a : zero(a)

function grade(a::MixedMultivector{sig}, k::Integer) where sig
	n = dimension(sig)
	i = bits_to_index(a, bits_first_of_grade(k))
	Multivector{sig,k}(a.components[i:i+binomial(n, k) - 1])
end

# experimental notation for even/odd grade projections
# grade(a::MixedMultivector, ::typeof(+)) = sum(grade(a, k) for k ∈ 0:2:dimension(a))
# grade(a::MixedMultivector, ::typeof(-)) = sum(grade(a, k) for k ∈ 1:2:dimension(a))

scalar(a::AbstractMultivector) = getcomponent(a)

isscalar(a::HomogeneousMultivector) = iszero(grade(a))
isscalar(a::MixedMultivector) = iszero(a.components[begin + 1:end])

pseudoscalar(a::AbstractMultivector) = getcomponent(a, (1:dimension(a))...)
pseudoscalar(a::MixedMultivector) = a.components[end]

ispseudoscalar(a::HomogeneousMultivector{sig,k}) where {sig,k} = k == dimension(sig)
ispseudoscalar(a::MixedMultivector) = iszero(a.components[begin:end - 1])


#= AUTOMORPHISMS & DUALITY OPERATIONS =#

reversionsign(k, x=1) = mod(k, 4) <= 1 ? +x : -x

"""
	~(a::AbstractMultivector)
	reversion(a::AbstractMultivector)

The reverse of a multivector.

The reversion antiautomorphism ``ρ``is the linear map defined by
``ρ(u) = u`` on vectors and ``ρ(ab) = ρ(b)ρ(a)`` on higher-grade elements.

See also [`involute`](@ref).
"""
reversion, ~

reversion(a::Scalar) = a
reversion(a::HomogeneousMultivector{sig,k}) where {sig,k} = reversionsign(k, a)
reversion(a::MixedMultivector) = mapcomponents(b -> reversion(b).coeff, a)

Base.:~(a::AbstractMultivector) = reversion(a)

"""
	involute(a::AbstractMultivector)

The involute of a multivector, negating the odd-grade part.

The involution automorphism ``ι`` is the linear map defined by
``ι(u) = -u`` on vectors and ``ι(ab) = ι(a)ι(b)`` on higher-grade elements.

See also [`reversion`](@ref).
"""
involute

involute(a::Scalar) = a
involute(a::HomogeneousMultivector) = iseven(grade(a)) ? a : -a
involute(a::MixedMultivector) = mapcomponents(u -> (iseven(grade(u)) ? u : -u).coeff, a)

# TODO: show the Clifford conjugate cong = reversion∘involute be defined as the adjoint?



unit_pseudoscalar(a::AbstractMultivector) = Blade{signature(a)}(one(eltype(a)), bits_first_of_grade(dimension(a)))

# TODO: which type of dual operation? aI, Ia, a/I, I\a?
# the Hodge dual is defined as ``φ ∧ ⋆ϕ = ⟨φ,ϕ⟩ vol`` for forms of the same degree
# so the choice `dual(a) = aI` agrees with the Hodge dual on blades of the same grade.
dual(a::AbstractMultivector) = a*unit_pseudoscalar(a)
dual(a::HomogeneousMultivector) = a⨼unit_pseudoscalar(a)


# experimental

flip_bits(a::Blade) = Blade{signature(a)}(a.coeff, bits_first_of_grade(dimension(a)) ⊻ a.bits)
flip_bits(a::Multivector) = best_type(a; grade=Val(dimension(a) - grade(a)))(reverse(a.components))
flip_bits(a::MixedMultivector) = best_type(a)(reverse(a.components))