#= Equality =#

Base.:(==)(a::Blade{Sig}, b::Blade{Sig}) where Sig = bitsof(a) == bitsof(b) ? a.coeff == b.coeff : iszero(a) && iszero(b)
Base.:(==)(a::Multivector{Sig}, b::Multivector{Sig}) where {Sig} = grade(a) == grade(b) ? all(a.components .== b.components) : iszero(a) && iszero(b)
Base.:(==)(a::MixedMultivector{Sig}, b::MixedMultivector{Sig}) where {Sig} = all(a.components .== b.components)

Base.:(==)(a::AbstractMultivector, b::Number) = iszero(b) && iszero(a) || isscalar(a) && scalarpart(a) == b
Base.:(==)(a::Number, b::AbstractMultivector) = iszero(a) && iszero(b) || isscalar(b) && a == scalarpart(b)

# equality between different multivector types
# TODO: implement without conversions
Base.:(==)(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig} = ==(largest_type(a, b).((a, b))...)


#= Scalar Multiplication =#

scalar_multiply(a::Blade, b) = Blade{signature(a)}(bitsof(a) => a.coeff*b)
scalar_multiply(a, b::Blade) = Blade{signature(b)}(bitsof(b) => a*b.coeff)

scalar_multiply(a::CompositeMultivector, b) = constructor(a)(a.components*b)
scalar_multiply(a, b::CompositeMultivector) = constructor(b)(a*b.components)

Base.:*(a::AbstractMultivector, b::Number) = scalar_multiply(a, b)
Base.:*(a::Number, b::AbstractMultivector) = scalar_multiply(a, b)
Base.:-(a::AbstractMultivector) = -one(eltype(a))*a

Base.:/(a::AbstractMultivector, b::Number) = a*inv(promote_to(eltype(a), b))
Base.:\(a::Number, b::AbstractMultivector) = inv(promote_to(eltype(b), a))*b

Base.://(a::AbstractMultivector, b::Number) = a*(one(b)//b)



#= Addition =#

add!(a::Multivector, b::Blade) = (a.components[mv_index(b)] += b.coeff; a)
add!(a::MixedMultivector, b::Blade) = (a.components[mmv_index(b)] += b.coeff; a)

add!(a::Multivector, b::Multivector) = (a.components += b.components; a)
add!(a::MixedMultivector, b::MixedMultivector) = (a.components += b.components; a)

function add!(a::MixedMultivector, b::Multivector)
	offset = multivector_index_offset(grade(b), dimension(b))
	a.components[mmv_slice(b)] = b.components
	a
end

# add alike types
Base.:+(As::Multivector{Sig,K}...) where {Sig,K} = Multivector{Sig,K}(sum(a.components for a ∈ As))
Base.:+(As::MixedMultivector{Sig}...) where {Sig} = MixedMultivector{Sig}(sum(a.components for a ∈ As))

# convert unalike to alike
Base.:+(As::HomogeneousMultivector{Sig,K}...) where {Sig,K} = +(Multivector.(As)...)
Base.:+(As::AbstractMultivector{Sig}...) where {Sig} = +(MixedMultivector.(As)...)

Base.:-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)




#= Geometric Multiplication =#

function geometric_prod(a::Blade{Sig}, b::Blade{Sig}) where {Sig}
	factor, bits = geometric_prod_bits(Sig, bitsof(a), bitsof(b))
	Blade{Sig}(bits => factor*(a.coeff*b.coeff))
end

function geometric_prod(a::AbstractMultivector{Sig}, b::AbstractMultivector{Sig}) where {Sig}
	T = promote_type(eltype(a), eltype(b))
	C = with_eltype(componentstype(Sig), T)
	ab = zero(MixedMultivector{Sig,C})
	for (abits, acoeff) ∈ nonzero_components(a), (bbits, bcoeff) ∈ nonzero_components(b)
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		i = bits_to_mmv_index(bits, dimension(Sig))
		ab.components[i] += factor*(acoeff*bcoeff)
	end
	ab
end

Base.:*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)