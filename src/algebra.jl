#= Equality =#

Base.:(==)(a::Blade{Sig}, b::Blade{Sig}) where Sig = bitsof(a) == bitsof(b) ? a.coeff == b.coeff : iszero(a) && iszero(b)
Base.:(==)(a::(Multivector{Sig,K,S} where K), b::(Multivector{Sig,K,S} where K)) where {Sig,S} = grade(a) == grade(b) ? a.components == b.components : iszero(a) && iszero(b)
Base.:(==)(a::MixedMultivector{Sig,S}, b::MixedMultivector{Sig,S}) where {Sig,S} = a.components == b.components

Base.:(==)(a::AbstractMultivector, b::Number) = iszero(b) && iszero(a) || isscalar(a) && scalar(a) == b
Base.:(==)(a::Number, b::AbstractMultivector) = iszero(a) && iszero(b) || isscalar(b) && a == scalar(b)



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