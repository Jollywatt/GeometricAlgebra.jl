#= SCALAR MULTIPLICATION =#

*(a::Blade{sig,k,bits}, b::Scalar) where {sig,k,bits} = Blade{sig,k,bits}(a.coeff*b)
*(a::Scalar, b::Blade{sig,k,bits}) where {sig,k,bits} = Blade{sig,k,bits}(a*b.coeff)

constructor(::Multivector{sig,k}) where {sig,k} = Multivector{sig,k}
constructor(::MixedMultivector{sig}) where sig = MixedMultivector{sig}
*(a::CompositeMultivector{<:AbstractVector}, b::Scalar) = constructor(a)(a.components*b)
*(a::Scalar, b::CompositeMultivector{<:AbstractVector}) = constructor(b)(a*b.components)

-(a::AbstractMultivector) = -one(eltype(a))*a

/(a::AbstractMultivector, b::Scalar) = a*inv(b)
\(a::Scalar, b::AbstractMultivector) = inv(a)*b



#= ADDITION =#

# mutating addition which does not require same element/storage types
function add!(a::CompositeMultivector, b::Blade)
	a.components[bits_to_key(a, bitsof(b))] += b.coeff
	a
end
function add!(a::CompositeMultivector, b::CompositeMultivector)
	for u ∈ blades(b)
		add!(a, u)
	end
	a
end

function add(As::(HomogeneousMultivector{sig,k} where sig)...) where k
	Σ = zero(best_type(Multivector, As..., grade=Val(k)))
	for a ∈ As, b ∈ blades(a)
		add!(Σ, b)
	end
	Σ
end
function add(As::AbstractMultivector...)
	Σ = zero(best_type(MixedMultivector, As...))
	for a ∈ As, b ∈ blades(a)
		add!(Σ, b)
	end
	Σ
end

add(a::Union{Blade,CompositeMultivector}) = a

+(a::AbstractMultivector...) = add(a...)
-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)

# scalar addition
+(a::AbstractMultivector{sig}, b::Scalar) where sig = a + Blade{sig}(b, bits_scalar())
+(a::Scalar, b::AbstractMultivector) = b + a
-(a::AbstractMultivector, b::Scalar) = a + (-b)
-(a::Scalar, b::AbstractMultivector) = a + (-b)




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
	Blade{sig,grade(bits),bits}(factor*(a.coeff*b.coeff))
end

function geometric_prod(a::AbstractMultivector, b::AbstractMultivector)
	Π = zero(best_type(MixedMultivector, a, b))
	for u ∈ blades(a), v ∈ blades(b)
		add!(Π, geometric_prod(u, v))
	end
	Π
end

"""
	*(::AbstractMultivector, ::AbstractMultivector)

Geometric product. See `$(fullname(geometric_prod))`.
"""
*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)



#= GRADED PRODUCTS =#

"""
	homogeneous_prod(a, b, Val(k))

The grade `k` part of the geometric product of `a` and `b`.
Equivalent to, but more efficient than, `grade(geometric_prod(a, b), k)`.
"""
homogeneous_prod

function homogeneous_prod(a::Blade, b::Blade, ::Val{k}) where k
	bits = bitsof(a) ⊻ bitsof(b)
	if grade(bits) == k
		geometric_prod(a, b)
	else
		zero(best_type(Blade, a, b; bits=Val(bits)))
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
			add!(Π, homogeneous_prod(u, v, k))
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
		add!(Π, homogeneous_prod(u, v, k))
	end
	Π
end
function graded_prod(a::HomogeneousMultivector, b::HomogeneousMultivector, grade_selector::Function)
	homogeneous_prod(a, b, grade_selector(grade(a), grade(b)))
end

∧(a::AbstractMultivector, b::AbstractMultivector) = graded_prod(a, b, +)
⋅(a::AbstractMultivector, b::AbstractMultivector) = graded_prod(a, b, abs∘-)



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

# returns (2^dim)-element vector of components for every basis blade in the algebra
full_components_vector(a::MixedMultivector{sig,<:AbstractVector}) where sig = a.components
function full_components_vector(a::HomogeneousMultivector)
	fcv = zeros(eltype(a), 2^dimension(a)) # TODO: use sparse vector? fcv is large and will remain mostly empty
	for b ∈ blades(a)
		fcv[begin + bitsof(b)] = b.coeff
	end
	fcv
end

function inv_matrixmethod(a::AbstractMultivector)
	A = hcat([full_components_vector(a*blade_like(a, 1, unsigned(i - 1)))
		for i ∈ 1:2^dimension(a)]...)
	id = full_components_vector(one(a))
	inv_components = A\id
	a⁻¹ = MixedMultivector{signature(a)}(inv_components)
	# result should be converted into the same storage type as original
end


/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b

/(a::Scalar, b::AbstractMultivector) = a*inv(b)
\(a::AbstractMultivector, b::Scalar) = inv(a)*b



#= EXPONENTIATION =#

# if a² is a scalar, then a²ⁿ = |a²|ⁿ, a²ⁿ⁺¹ = |a²|ⁿa
function power_with_scalar_square(a::AbstractMultivector, a², p::Integer)
	# if p is even, p = 2n; if odd, p = 2n + 1
	aⁿ = a²^fld(p, 2)
	iseven(p) ? aⁿ*one(a) : aⁿ*a
end

function power_by_squaring(a::AbstractMultivector, p::Integer)
	p >= 0 || return power_by_squaring(a, abs(p))
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
	a² = a*a
	if isscalar(a²)
		power_with_scalar_square(a, scalar(a²), p)
	else
		power_by_squaring(a, p)
	end
end



#= GRADE PROJECTION =#

grade(a::Scalar) = 0
grade(a, k::Integer) where sig = grade(a) == k ? a : zero(a)
function grade(a::MixedMultivector{sig,S}, k::Integer) where {sig,S}
	b = zero(Multivector{sig,k,S})
	for u ∈ blades(a)
		if grade(u) == k
			add!(b, u)
		end
	end
	b
end

# experimental notation for even/odd grade projections
# grade(a::MixedMultivector, ::typeof(+)) = sum(grade(a, k) for k ∈ 0:2:dimension(a))
# grade(a::MixedMultivector, ::typeof(-)) = sum(grade(a, k) for k ∈ 1:2:dimension(a))

grades(a::MixedMultivector) = sort(unique(grade.(Iterators.filter(!iszero, blades(a)))))

scalar(a::AbstractMultivector) = getcomponent(a)

isscalar(a::HomogeneousMultivector) = iszero(grade(a))
isscalar(a::MixedMultivector) = all(iszero.(grades(a)))



#= AUTOMORPHISMS & DUALITY OPERATIONS =#

reversionsign(k, x=1) = mod(k, 4) <= 1 ? +x : -x

"""
	~(a::AbstractMultivector)
	reversion(a::AbstractMultivector)

The reverse of a multivector.

The reversion antiautomorphism ``ρ`` is defined by ``ρ(u) = u`` on vectors
and ``ρ(ab) = ρ(b)ρ(a)`` on higher-grade elements.
"""
reversion, ~

reversion(a::Scalar) = a
reversion(a::HomogeneousMultivector{sig,k}) where {sig,k} = reversionsign(k, a)
reversion(a::MixedMultivector) = mapcomponents(b -> reversion(b).coeff, a)

Base.:~(a::AbstractMultivector) = reversion(a)

"""
	involute(a::AbstractMultivector)

The involute of a multivector, negating the odd-grade part.

The involution automorphism ``ι`` is defined by ``ι(u) = -u`` on vectors
and ``ι(ab) = ι(a)ι(b)`` on higher-grade elements.
"""
involute

involute(a::Scalar) = a
involute(a::HomogeneousMultivector) = iseven(grade(a)) ? a : -a
involute(a::MixedMultivector) = mapcomponents(u -> (iseven(grade(u)) ? u : -u).coeff, a)

# TODO: show the Clifford conjugate cong = reversion∘involute be defined as the adjoint?

