
#= SCALAR MULTIPLICATION =#


*(a::Blade{sig,k,bits}, b::Scalar) where {sig,k,bits} = Blade{sig,k,bits}(a.coeff*b)
*(a::Scalar, b::Blade{sig,k,bits}) where {sig,k,bits} = Blade{sig,k,bits}(a*b.coeff)

constructor(::Multivector{sig,k}) where {sig,k} = Multivector{sig,k}
constructor(::MixedMultivector{sig}) where sig = MixedMultivector{sig}
*(a::CompositeMultivector{<:AbstractVector}, b::Scalar) = constructor(a)(a.components*b)
*(a::Scalar, b::CompositeMultivector{<:AbstractVector}) = constructor(b)(a*b.components)



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
function homogeneous_prod(a::HomogeneousMultivector{k1}, b::HomogeneousMultivector{k2}, ::Val{k}) where {k1,k2,k}
	Π = zero(best_type(Multivector, a, b; grade=Val(k)))
	k ∈ abs(k1 - k2):2:(k1 + k2) && return Π
	for u ∈ blades(a), v ∈ blades(b)
		add!(Π, homogeneous_prod(u, v, k))
	end
	Π
end
function homogeneous_prod(a::AbstractMultivector, b::AbstractMultivector, ::Val{k}) where k
	Π = zero(best_type(Multivector, a, b; grade=Val(k)))
	0 <= k <= dimension(a) || return Π
	for u ∈ blades(a), v ∈ blades(b)
		add!(Π, homogeneous_prod(u, v, k))
	end
	Π
end
homogeneous_prod(a, b, k::Integer) = homogeneous_prod(a, b, Val(k))

"""
	graded_prod(a, b, grade_selector::Function)

Homogenised geometric product with resulting grade given as a function of the grades of `a` and `b`.

If ``a_p`` and ``b_q`` are multivectors of grade ``p`` and ``q`` respectively, then this is
equal to ``⟨a_p b_q⟩_{k(p,q)}`` where ``k(p,q)`` is given by `grade_selector`.
If `a` or `b` are multi-graded, then this extends linearly to ``∑_{pq} ⟨⟨a⟩_p ⟨b⟩_q⟩_{k(p,q)}``.
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





#= COMPONENT ACCESS =#

bits_to_key(a::Multivector{sig,k,<:AbstractVector}, bits::Unsigned) where {sig,k} = bits_to_linear_index(bits)
bits_to_key(a::MixedMultivector{sig,<:AbstractVector}, bits::Unsigned) where sig = firstindex(a.components) + bits
bits_to_key(a::CompositeMultivector{<:AbstractDict}, bits::Unsigned) = bits

parity_sign(I) = iseven(parity(sortperm(collect(I)))) ? +1 : -1

function getcomponent(a::Blade, I...)
	indices_to_bits(I) == bitsof(a) || return zero(eltype(a))
	s = parity_sign(I)
	s*a.coeff
end

getcomponent(a::Multivector, I...) = zero(eltype(a))
function getcomponent(a::Multivector{sig,k}, I::Vararg{Any,k}) where {sig,k}
	s = parity_sign(I)
	s*a.components[bits_to_key(a, indices_to_bits(I))]
end

function getcomponent(a::MixedMultivector, I...)
	s = parity_sign(I)
	s*a.components[bits_to_key(a, indices_to_bits(I))]
end

Base.getindex(a::AbstractMultivector, I::Integer...) = getcomponent(a, I...)
Base.getindex(a::AbstractMultivector) = getcomponent(a)

# experimental notations for grade selection
Base.getindex(a::AbstractMultivector, I::Vector{<:Integer}) = grade(a, Iterators.only(I))
# Base.getindex(a::AbstractMultivector; grade) = GeometricAlgebra2.grade(a, grade)


function mapcomponents(f, a::CompositeMultivector)
	a′ = zero(a)
	for b ∈ blades(a)
		a′.components[bits_to_key(a, bitsof(b))] = f(b)
	end
	a′
end



#= GRADE PROJECTION =#

grade(a::Blade{sig}, k::Integer) where sig = grade(a) == k ? a : zero(a)
grade(a::Multivector{sig}, k::Integer) where sig = grade(a) == k ? a : zero(a)
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
grade(a::MixedMultivector, ::typeof(+)) = sum(grade(a, k) for k ∈ 0:2:dimension(a))
grade(a::MixedMultivector, ::typeof(-)) = sum(grade(a, k) for k ∈ 1:2:dimension(a))

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
involute(a::HomogeneousMultivector) = iseven(grade(a)) ? a : -a
involute(a::MixedMultivector) = mapcomps(u -> (iseven(grade(u)) ? u : -u).coeff, a)

# TODO: show the Clifford conjugate cong = reversion∘involute be defined as the adjoint?



-(a::AbstractMultivector) = -one(eltype(a))*a
-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)

+(a::AbstractMultivector{sig}, b::Scalar) where sig = a + Blade{sig}(b, bits_scalar())
+(a::Scalar, b::AbstractMultivector) = b + a
-(a::AbstractMultivector, b::Scalar) = a + (-b)
-(a::Scalar, b::AbstractMultivector) = a + (-b)