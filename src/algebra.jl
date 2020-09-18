# SCALAR MULTIPLICATION
# must respect the possibility of non-commutativity!

compmul(comps::AbstractVector, x::Number) = comps.*x
compmul(x::Number, comps::AbstractVector) = x.*comps
compmul(comps::AbstractDict, x::Number) = Dict(ublade => a*val for (ublade, val) ∈ b.comps)
compmul(x::Number, comps::AbstractDict) = Dict(ublade => val*b for (ublade, val) ∈ a.comps)

*(a::Number, b::Blade) = Blade{signature(b)}(a*b.coeff, b.ublade)
*(a::Blade, b::Number) = Blade{signature(a)}(a.coeff*b, a.ublade)

*(a::Number, b::Multivector) = Multivector{signature(b)}(grade(b), compmul(a, b.comps))
*(a::Multivector, b::Number) = Multivector{signature(a)}(grade(a), compmul(a.comps, b))

*(a::Number, b::MixedMultivector) = MixedMultivector{signature(b)}(compmul(a, b.comps))
*(a::MixedMultivector, b::Number) = MixedMultivector{signature(a)}(compmul(a.comps, b))

-(a::Blade) = typeof(a)(-a.coeff, a.ublade)
-(a::AbstractMultivector) = -one(eltype(a))*a


# (MULTI)VECTOR ADDITION

function +(as::HomogeneousMultivector{k}...) where k
	T = best_type(Multivector, as...; k)
	Σ = zero(T)
	for a ∈ as, (ublade, val) ∈ comps(a)
		addcomp!(Σ, ublade, val)
	end
	Σ
end

function +(as::AbstractMultivector...)
	T = best_type(MixedMultivector, as...)
	Σ = zero(T)
	for a ∈ as, (ublade, val) ∈ comps(a)
		addcomp!(Σ, ublade, val)
	end
	Σ
end


function +(a::MixedMultivector, b::Number)
	ab = deepcopy(a)
	ab[] += b
	ab
end
+(a::Number, b::MixedMultivector) = +(b, a)


+(a::AbstractMultivector, b::Number) = convert(MixedMultivector, a) + b
+(a::Number, b::AbstractMultivector) = a + convert(MixedMultivector, b)


-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)
-(a::AbstractMultivector, b::Number) = a + (-b)
-(a::Number, b::AbstractMultivector) = a + (-b)





#= ALGEBRAIC PRODUCTS =#

function geometric_prod(a1::AbstractMultivector, a2::AbstractMultivector)
	T = best_type(MixedMultivector, a1, a2)
	sig = signature(T)
	a12 = zero(T)
	for (u1, v1) ∈ comps(a1), (u2, v2) ∈ comps(a2)
		factor, u12 = ubladeprod(sig, u1, u2)
		v12 = v1*v2
		addcomp!(a12, u12, factor*v12)
	end
	a12
end
function homogeneous_prod(a1::AbstractMultivector, a2::AbstractMultivector, k)
	T = best_type(Multivector, a1, a2; k)
	sig = signature(T)
	a12 = zero(T)
	for (u1, v1) ∈ comps(a1), (u2, v2) ∈ comps(a2)
		u12 = ublade_xor(u1, u2)
		if ublade_grade(u12) == k
			factor, _ = ubladeprod(sig, u1, u2)
			v12 = v1*v2
			addcomp!(a12, u12, factor*v12)
		end
	end
	a12
end
function graded_prod(a1::AbstractMultivector, a2::AbstractMultivector, grade_selector)
	T = best_type(MixedMultivector, a1, a2)
	sig = signature(T)
	a12 = zero(T)
	for (u1, v1) ∈ comps(a1), (u2, v2) ∈ comps(a2)
		u12 = ublade_xor(u1, u2)
		if grade_selector(ublade_grade(u1), ublade_grade(u2)) == ublade_grade(u12)
			factor, _ = ubladeprod(sig, u1, u2)
			v12 = v1*v2
			addcomp!(a12, u12, factor*v12)
		end
	end
	a12
end

# specialisations for blades
function geometric_prod(a1::Blade, a2::Blade)
	sig = shared_sig(a1, a2)
	factor, u12 = ubladeprod(sig, a1.ublade, a2.ublade)
	v12 = a1.coeff*a2.coeff
	Blade{sig}(factor*v12, u12)
end
homogeneous_prod(a::Blade, b::Blade, k) = grade(geometric_prod(a, b), k)

graded_prod(a::HomogeneousMultivector, b::HomogeneousMultivector, grade_selector) = homogeneous_prod(a, b, grade_selector(grade(a), grade(b)))

*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)

"""
```
a⋅b
dot(a, b)
```

"Fat" dot product of multivectors, giving the lowest-grade part of the geometric product.
Defined by ``a⋅b = ∑_{ij} ⟨⟨a⟩_i * ⟨b⟩_j⟩_|i - j|``.
"""
dot(a, b) = graded_prod(a, b, abs∘-)
const ⋅ = dot

"""
```
a∧b
wedge(a, b)
```

Wedge product of multivectors, giving the highest-grade part of the geometric product `a*b`.
Defined by ``a∧b = ∑_{ij} ⟨⟨a⟩_i * ⟨b⟩_j⟩_(i + j)``.
"""
wedge(a, b) = graded_prod(a, b, +)
const ∧ = wedge

"""
```
a∗b
scalar_prod(a, b)
```

Scalar product of multivectors, equivalent to `grade(a*b, 0)`.
"""
scalar_prod(a, b) = homogeneous_prod(a, b, 0)
const ∗ = scalar_prod

contractr(a, b) = graded_prod(a, b, (r, s) -> r - s)
const ⨽ = contractr

contractl(a, b) = graded_prod(a, b, (r, s) -> s - r)
const ⨼ = contractl

"""
```
a×b
commutator_prod(a, b)
```

Commutator product of multivectors, given by ``a × b = (a*b - b*a)/2``.

"""
commutator_prod(a, b) = (a*b - b*a)/2




# REVERSION, INVOLUTION, DUALITY, ETC

# @pure ?
reversionsign(k::Integer, c) = mod(k, 4) <= 1 ? c : -c
reversionsign(k::Integer) = reversionsign(k, 1)

reversion(a::Number) = a
reversion(a::HomogeneousMultivector{k}) where k = reversionsign(k, a)
reversion(a::MixedMultivector) = mapcomp((u, v) -> reversionsign(ublade_grade(u), v), a)
~(a::AbstractMultivector) = reversion(a)

involute(a::HomogeneousMultivector) = iseven(grade(a)) ? a : -a
involute(a::MixedMultivector) = mapcomp((u, v) -> iseven(ublade_grade(u)) ? v : -v, a)


"""
```
★(a)
hodgedual(a)
```

Hodge dual or 'star' operator, satisfying
``a∧★(b) == (a∗b)vol(a)``.(?)
"""
hodgedual(a::AbstractMultivector) = reversion(a)*vol(a)
const ★ = hodgedual



# MULTIPLICATIVE INVERSES

Base.inv(a::Blade) = a/(a*a)[] # blades are guaranteed to have scalar squares
function Base.inv(a::AbstractMultivector)
	ã = reversion(a)
	square = ã*a
	isscalar(square) || error("multivector is not invertible")
	ã/square[]
end



/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
/(a::AbstractMultivector, b::Number) = a*inv(b)
/(a::Number, b::AbstractMultivector) = a*inv(b)

\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b
\(a::AbstractMultivector, b::Number) = inv(a)*b
\(a::Number, b::AbstractMultivector) = inv(a)*b


# EXPONENTIATION

function ^(a::Blade, p::Integer)
	p >= 0 || return inv(a)^abs(p)
	blade = iseven(p) ? one(a) : oneunit(a)
	square, _ = ubladeprod(signature(a), a.ublade, a.ublade)
	# square = 3
	factor = square^fld(p, 2) # floored division
	coeff = factor*a.coeff^p
	coeff*blade
end

function powbysquaring(a::AbstractMultivector, p::Integer)
	p >= 0 || return powbysquaring(inv(m), abs(p))
	Π = one(best_type(MixedMultivector, a))
	aᵖ = a
	while p > 0
		if !iszero(p & 1)
			Π *= aᵖ
		end
		aᵖ *= aᵖ
		p >>= 1
	end
	Π
end

^(a::AbstractMultivector, p::Integer) = powbysquaring(a, p)
# ^(a::AbstractMultivector, p) = error("unsupported multivector exponent type $(typeof(p))")

function Base.exp(a::Blade)
	factor, _ = ubladeprod(signature(a), a.ublade, a.ublade)
	coeff = a.coeff
	s = sign(factor)
	if iszero(s)
		one(a) + a
	elseif s > 0
		cosh(coeff) + oneunit(a)*sinh(coeff)
	else
		cos(coeff) + oneunit(a)*sin(coeff)
	end
end


_sumcoeff(a::Blade) = a.coeff
_sumcoeff(a::CompositeMultivector{sig,<:AbstractVector}) where sig = sum(a.comps)
_sumcoeff(a::CompositeMultivector{sig,<:AbstractDict}) where sig = sum(values(a.comps))

# what's the proper way to implement a stable exp which works with BigFloats?
function Base.exp(a::AbstractMultivector, n=20)
	eᵃ = one(best_type(MixedMultivector, a))
	term = one(eᵃ)
	for i ∈ 1:n
		term *= a/eltype(a)(i)
		eᵃ += term
	end
	eᵃ
end