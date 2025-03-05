struct ProductNode{K}
	x::WeightedSet{K,Int}
end

struct SumNode{K,V}
	x::WeightedSet{ProductNode{K},V}
end

const Π = ProductNode
const Σ = SumNode

Base.copy(a::Union{Π,Σ}) = typeof(a)(copy(a.x))

(a::Π == b::Π) = a.x == b.x
(a::Σ == b::Σ) = a.x == b.x
Base.hash(a::Union{Π,Σ}, seed::UInt) = hash(a.x, seed)

Π(ps::Pair...) = Π(WeightedSet(ps...))
Σ(ps::Pair...) = Σ(WeightedSet(ps...))
Π(iter::Base.Generator) = Π(WeightedSet(iter))
Σ(iter::Base.Generator) = Σ(WeightedSet(iter))

Π{K}(iter::Base.Generator) where K = Π{K}(WeightedSet{K,Int}(iter))

Base.convert(::Type{Π{K}}, a::Π{K′}) where {K,K′} = Π(WeightedSet{K,Int}(a.x))
Σ(d::WeightedSet{Π,V}) where V = Σ(WeightedSet{Π{Any},V}(d))

Base.one(::Type{Π{K}}) where K = Π{K}(Dict{K,Int}())
Base.one(::Type{Σ{K,V}}) where {K,V} = Σ(one(Π{K}) => one(V))
Base.zero(::Type{Σ{K,V}}) where {K,V} = Σ{K,V}(Dict{Π{K},V}())

Base.one(a::Union{Π,Σ}) = one(typeof(a))
Base.zero(a::Union{Π,Σ}) = zero(typeof(a))

ismonomial(a::Π) = true
ismonomial(a::Σ) = isone(length(a.x))

#= multiplication of products =#

mul!(a::Π, b::Π) = Π(combine!(a.x, b.x))
(a::Π * b::Π) = Π(combine(a.x, b.x))

(a::Π / b::Π) = mul!(inv(b), a)
(a::Π \ b::Π) = mul!(inv(a), b)
(a::Union{Π,Σ} / b::Union{Π,Σ}) = a*inv(b)
(a::Union{Π,Σ} \ b::Union{Π,Σ}) = inv(a)*b

(a::Union{Π,Σ} * λ::Number) = scalarmul(a, λ)
(a::Union{Π,Σ} / λ::Number) = scalarmul(a, inv(λ))
(λ::Number * a::Union{Π,Σ}) = scalarmul(a, λ)
(λ::Number \ a::Union{Π,Σ}) = scalarmul(a, inv(λ))

(λ::Number / a::Π) = scalarmul(inv(a), λ)
(λ::Number / a::Σ) = scalarmul!(inv(a), λ)
(a::Π \ λ::Number) = scalarmul(inv(a), λ)
(a::Σ \ λ::Number) = scalarmul!(inv(a), λ)

Base.inv(a::Π) = Π(negate!(copy(a.x)))
(a::Π{K} ^ p::Integer) where K = Π{K}(k => p*v for (k, v) in a.x)


#= scalar multiplication of sums =#
scalarmul(a::Σ, λ::Number) = Σ(k => λ*v for (k, v) in a.x)
scalarmul!(a::Σ, λ::Number) = (mul!(a.x, λ); a)
-(a::Σ) = Σ(negate!(copy(a.x)))

#= scalar multiplication of products =#
scalarmul(a::Π, λ::Number) = Σ(a => λ)
+(a::Π) = Σ(a => +1)
-(a::Π) = Σ(a => -1)

#= addition of sums and products =#

add!(a::Σ, b::Σ) = (combine!(a.x, b.x); a)
add!(a::Σ, b::Π, λ) = (a.x[b] += λ; a)


(a::Π + b::Π) = Σ(a => 1, b => 1)
(a::Σ + b::Σ) = Σ(combine(a.x, b.x))
(a::Union{Π,Σ} + λ::Number) = a + λ*one(a)
(λ::Number + a::Union{Π,Σ}) = a + λ

function (a::Σ{K,V} + b::Π{K′}) where {K,V,K′}
	K >: K′ && return add!(copy(a), b, 1)
	c = zero(Σ{promote_type(K, K′),V})
	add!(c, a)
	add!(c, b, 1)
end
(a::Π + b::Σ) = b + a

(a::Union{Π,Σ} - b::Σ) = a + (-b)
(a::Σ - b::Π) = add!(copy(a), b, -1)
(a::Π - b::Π) = Σ(a => 1, b=> -1)
(λ::Number - a::Union{Π,Σ}) = (-a) + λ
(a::Union{Π,Σ} - λ::Number) = a + (-λ)


#= multiplication of sums =#

(a::Σ * b::Σ) = Σ(p*q => α*β for (p, α) in a.x, (q, β) in b.x)
(a::Σ * b::Π) = Σ(p*b => α for (p, α) in a.x)
(a::Π * b::Σ) = b*a


#= factor =#

function Base.gcd(a::Π, rest::Π...)
	common = copy(a)
	for b in rest
		for (k, v) in common.x
			common.x[k] = min(v, b.x[k])
		end
	end
	common
end

function firstcoeff(a::Σ)
	c = 0.0
	h = UInt(0)
	for (k, v) in a.x
		h′ = hash(k)
		h′ > h || continue
		h = h′
		c = float(v)
	end
	c
end

function fixmul!!(a::Σ{K,<:Integer}) where K
	c = firstcoeff(a)
	c, a/c
end
function fixmul!!(a::Σ{K,<:Real}) where K
	c = firstcoeff(a)
	c, scalarmul!(a, inv(c))
end

"""
	factor(x::SumNode)

Naively collect factors that are common to all terms in `x`.
For example, `x*y + x*z` becomes `x*(y + z)`, but `x^2 + 2x + 1` is left as is.
"""
function factor(a::Σ)
	iszero(a) && return a
	ismonomial(a) && return a
	commonfactor = reduce(gcd, keys(a.x))

	if isone(commonfactor)
		s, a = fixmul!!(a)
		Σ(Π(a => 1) => s/1)
	else
		s, rest = fixmul!!(a/commonfactor)
		Σ(Π(rest => 1)*commonfactor => s)
	end
end


#= exponentiation =#

function (a::Σ ^ p::Integer)
	iszero(a) && return a
	iszero(p) && return one(a)
	isone(p) && return a
	if ismonomial(a)
		k, v = only(a.x)
		Σ(k^p => v^p)
	else
		factor(a)^p
	end
end

# should return sum with float coeffs
function Base.inv(a::Σ)
	if ismonomial(a)
		k, v = only(a.x)
		Σ(inv(k) => inv(v))
	else
		inv(factor(a))
	end
end


#= special functions =#

for fn in [:sqrt, :cosh, :sinh, :cos, :sin]
	@eval Base.$fn(a::Union{Π,Σ}) = Π(Expr(:call, $(Meta.quot(fn)), toexpr(a)) => 1)
end

#= interaction with vectors =#

(a::Union{Π,Σ} * v::AbstractArray{<:Union{Π,Σ}}) = Ref(a).*v
(v::AbstractArray{<:Union{Π,Σ}} * a::Union{Π,Σ}) = v.*Ref(a)
