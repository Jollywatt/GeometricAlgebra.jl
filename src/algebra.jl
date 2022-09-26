#= Equality =#

Base.:(==)(a::Blade{Sig}, b::Blade{Sig}) where Sig = bitsof(a) == bitsof(b) ? a.coeff == b.coeff : iszero(a) && iszero(b)
Base.:(==)(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig} = grade(a) == grade(b) ? a.comps == b.comps : iszero(a) && iszero(b)
Base.:(==)(a::MixedMultivector{Sig}, b::MixedMultivector{Sig}) where {Sig} = a.comps == b.comps

Base.:(==)(a::AbstractMultivector, b::Number) = isscalar(a) && scalarpart(a) == b
Base.:(==)(a::Number, b::AbstractMultivector) = isscalar(b) && a == scalarpart(b)

# equality between different multivector types
# TODO: implement without conversions?
Base.:(==)(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig} = let T = largest_type(a, b)
	T(a) == T(b)
end



#= Approximate Equality =#

isapproxzero(a; kwargs...) = isapprox(a, zero(a); kwargs...)
isapproxzero(a::Blade; kwargs...) = isapproxzero(a.coeff; kwargs...)
isapproxzero(a::CompositeMultivector; kwargs...) = isapproxzero(a.comps; kwargs...)

Base.isapprox(a::Blade{Sig}, b::Blade{Sig}; kwargs...) where Sig = bitsof(a) == bitsof(b) ? isapprox(a.coeff, b.coeff; kwargs...) : isapproxzero(a) && isapproxzero(b)
Base.isapprox(a::Multivector{Sig}, b::Multivector{Sig}; kwargs...) where {Sig} = grade(a) == grade(b) ? isapprox(a.comps, b.comps; kwargs...) : isapproxzero(a) && isapproxzero(b)
Base.isapprox(a::MixedMultivector{Sig}, b::MixedMultivector{Sig}; kwargs...) where {Sig} = isapprox(a.comps, b.comps; kwargs...)

# promote scalar to target multivector type and compare component arrays
Base.:isapprox(a::Blade, b::Number; kwargs...) = isapprox(Multivector(a), b; kwargs...)
Base.:isapprox(a::CompositeMultivector, b::Number; kwargs...) = isapprox(a, zero(a) + b; kwargs...)
Base.:isapprox(a::Number, b::AbstractMultivector; kwargs...) = isapprox(b, a; kwargs...)

Base.isapprox(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}; kwargs...) where {Sig} = let T = largest_type(a, b)
	isapprox(T(a), T(b); kwargs...)
end



#= Scalar Multiplication =#

const Scalar = Union{Number,SymbolicUtils.Symbolic}

scalar_multiply(a::Blade{Sig,K}, b) where {Sig,K} = Blade{Sig,K}(bitsof(a) => a.coeff*b)
scalar_multiply(a, b::Blade{Sig,K}) where {Sig,K} = Blade{Sig,K}(bitsof(b) => a*b.coeff)

scalar_multiply(a::CompositeMultivector, b) = constructor(a)(a.comps*b)
scalar_multiply(a, b::CompositeMultivector) = constructor(b)(a*b.comps)

Base.:*(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
Base.:*(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)
Base.:-(a::AbstractMultivector) = -realone(eltype(a))*a

promote_to(T, x) = convert(promote_type(T, typeof(x)), x)
Base.:/(a::AbstractMultivector, b::Scalar) = a*inv(promote_to(eltype(a), b))
Base.:\(a::Scalar, b::AbstractMultivector) = inv(promote_to(eltype(b), a))*b

Base.://(a::AbstractMultivector, b::Scalar) = a*(one(b)//b)



#= Addition =#

add!(a::Multivector, b::Blade) = (a.comps[mv_index(b)] += b.coeff; a)
add!(a::MixedMultivector, b::Blade) = (a.comps[mmv_index(b)] += b.coeff; a)

add!(a::Multivector, b::Multivector) = (a.comps[:] += b.comps; a)
add!(a::MixedMultivector, b::MixedMultivector) = (a.comps[:] += b.comps; a)

function add!(a::MixedMultivector, b::Multivector)
	offset = multivector_index_offset(grade(b), dimension(b))
	a.comps[mmv_slice(b)] = b.comps
	a
end

# add alike types
Base.:+(As::Multivector{Sig,K}...) where {Sig,K} = Multivector{Sig,K}(sum(a.comps for a ∈ As))
Base.:+(As::MixedMultivector{Sig}...) where {Sig} = MixedMultivector{Sig}(sum(a.comps for a ∈ As))

# convert unalike to alike # TODO: reduce intermediate allocations
Base.:+(As::HomogeneousMultivector{Sig,K}...) where {Sig,K} = +(Multivector.(As)...)
Base.:+(As::AbstractMultivector{Sig}...) where {Sig} = +(MixedMultivector.(As)...)

Base.:-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)



#= Scalar Addition =#

add_scalar(a::AbstractMultivector{Sig}, b) where {Sig} = a + Blade{Sig}(0 => b)
function add_scalar(a::MixedMultivector, b)
	constructor(a)(copy_setindex(a.comps, a.comps[1] + b, 1))
end

Base.:+(a::AbstractMultivector, b::Scalar) = add_scalar(a, b)
Base.:+(a::Scalar, b::AbstractMultivector) = add_scalar(b, a)

Base.:-(a::AbstractMultivector, b::Scalar) = add_scalar(a, -b)
Base.:-(a::Scalar, b::AbstractMultivector) = add_scalar(-b, a)



#= Geometric Multiplication =#

function geometric_prod(a::Blade{Sig}, b::Blade{Sig}) where {Sig}
	factor, bits = geometric_prod_bits(Sig, bitsof(a), bitsof(b))
	Blade{Sig}(bits => factor*(a.coeff*b.coeff))
end

function geometric_prod(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	c = zero(similar(MixedMultivector{Sig}, a, b))
	for (abits, ai) ∈ nonzero_components(a), (bbits, bi) ∈ nonzero_components(b)
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		i = bits_to_mmv_index(bits, dimension(Sig))
		c = setindex!!(c, c.comps[i] + factor*(ai*bi), i)
	end
	c
end

@generated geometric_prod_gen(a, b) = generated_multivector_function(geometric_prod, a, b)

geometric_prod(a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
geometric_prod(a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)
geometric_prod(a::Scalar, b::Scalar) = a*b

Base.:*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)



#= Derived Products =#

scalar_prod(a::Blade{Sig,K}, b::Blade{Sig,K}) where {Sig,K} = bitsof(a) == bitsof(b) ? scalarpart(a*b) : realzero(promote_type(eltype(a), eltype(b)))
scalar_prod(a::Blade{Sig}, b::Blade{Sig}) where {Sig} = realzero(promote_type(eltype(a), eltype(b)))

function scalar_prod(a::Multivector{Sig,K}, b::Multivector{Sig,K}) where {Sig,K}
	Blade{Sig}(0 => sum(geometric_square_factor.(Ref(Sig), bitsof(a)) .* (a.comps .* b.comps)))
end
scalar_prod(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig} = realzero(promote_type(eltype(a), eltype(b)))


function scalar_prod(a::MixedMultivector{Sig}, b::MixedMultivector{Sig}) where {Sig}
	Blade{Sig}(0 => sum(geometric_square_factor.(Ref(Sig), bitsof(a)) .* (a.comps .* b.comps)))
end
scalar_prod(a::AbstractMultivector, b::AbstractMultivector) = let T = largest_type(a, b)
	scalar_prod(T(a), T(b))
end



function graded_prod(grade_selector::Function, a::Blade{Sig}, b::Blade{Sig}) where {Sig}
	if count_ones(bitsof(a) ⊻ bitsof(b)) == grade_selector(grade(a), grade(b))
		a*b
	else
		Blade{Sig}(0 => zero(promote_type(eltype(a), eltype(b))))
	end
end

function graded_prod(grade_selector::Function, a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	c = zero(similar(MixedMultivector{Sig}, a, b))
	for (abits, ai) ∈ nonzero_components(a), (bbits, bi) ∈ nonzero_components(b)
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		if count_ones(bits) == grade_selector(count_ones(abits), count_ones(bbits))
			i = bits_to_mmv_index(bits, dimension(Sig))
			c = setindex!!(c, c.comps[i] + factor*(ai*bi), i)
		end
	end
	c
end

@generated function graded_prod_gen(grade_selector, a, b)
	generated_multivector_function(a, b) do a, b
		graded_prod(grade_selector.instance, a, b)
	end
end

# these are correct ONLY IF grade_selector(k, 0) == k == grade_selector(0, k)
graded_prod(grade_selector, a::AbstractMultivector, b::Scalar) = scalar_multiply(a, b)
graded_prod(grade_selector, a::Scalar, b::AbstractMultivector) = scalar_multiply(a, b)
graded_prod(grade_selector, a::Scalar, b::Scalar) = a*b

wedge(a, b) = graded_prod(+, a, b)
∧(a, b) = wedge(a, b)

inner(a, b) = graded_prod(abs∘-, a, b)
⋅(a, b) = inner(a, b)



#= Exponentiation =#

# if a² is a scalar, then a²ⁿ = |a²|ⁿ, a²ⁿ⁺¹ = |a²|ⁿa
function power_with_scalar_square(a, a², p::Integer)
	# if p is even, p = 2n; if odd, p = 2n + 1
	aⁿ = a²^fld(p, 2)
	iseven(p) ? aⁿ*one(a) : aⁿ*a
end

function power_by_squaring(a::CompositeMultivector{Sig,S}, p::Integer) where {Sig,S}
	p < 0 && return power_by_squaring(inv(a), abs(p))
	Π = one(MixedMultivector{Sig,S})
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

Base.:^(a::Blade, p::Integer) = power_with_scalar_square(a, scalarpart(a*a), p)
Base.literal_pow(::typeof(^), a::Blade{Sig}, ::Val{2}) where {Sig} = Blade{Sig,0}(0 => geometric_square_factor(Sig, bitsof(a))*a.coeff^2)

function Base.:^(a::CompositeMultivector{Sig,S}, p::Integer) where {Sig,S}
	# TODO: type stability?
	p == 0 && return one(a)
	p == 1 && return a
	a² = a*a
	p == 2 && return a²
	if isscalar(a²)
		power_with_scalar_square(a, scalarpart(a²), p)
	else
		power_by_squaring(a, p)
	end
end



#= Reversion =#

graded_multiply(f, a::Scalar) = f(0)*a
graded_multiply(f, a::HomogeneousMultivector) = f(grade(a))*a
function graded_multiply(f, a::MixedMultivector{Sig}) where Sig
	comps = copy(a.comps)
	dim = dimension(Sig)
	for k ∈ 0:dim
		comps[mmv_slice(Val(dim), Val(k))] *= f(k)
	end
	MixedMultivector{Sig}(comps)
end


reversion(a) = graded_multiply(reversion_sign, a)
Base.:~(a::AbstractMultivector) = reversion(a)


involution(a) = graded_multiply(k -> (-1)^k, a)

clifford_conj(a) = graded_multiply(a) do k
	(-1)^k*reversion_sign(k)
end