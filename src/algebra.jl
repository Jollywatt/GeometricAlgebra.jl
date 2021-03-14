#= SCALAR MULTIPLICATION =#
# must respect the possibility of non-commutativity!

compmul(comps::AbstractVector, x::Scalar) = comps.*Ref(x)
compmul(x::Scalar, comps::AbstractVector) = Ref(x).*comps
compmul(comps::AbstractDict, x::Scalar) = Dict(ublade => a*val for (ublade, val) ∈ b.comps)
compmul(x::Scalar, comps::AbstractDict) = Dict(ublade => val*b for (ublade, val) ∈ a.comps)

*(a::Scalar, b::Blade{sig,k}) where {sig,k} = Blade{sig,k}(a*b.coeff, b.ublade)
*(a::Blade{sig,k}, b::Scalar) where {sig,k} = Blade{sig,k}(a.coeff*b, a.ublade)

*(a::Scalar, b::Multivector{sig,k}) where {sig,k} = Multivector{sig,k}(compmul(a, b.comps))
*(a::Multivector{sig,k}, b::Scalar) where {sig,k} = Multivector{sig,k}(compmul(a.comps, b))

*(a::Scalar, b::MixedMultivector) = MixedMultivector{signature(b)}(compmul(a, b.comps))
*(a::MixedMultivector, b::Scalar) = MixedMultivector{signature(a)}(compmul(a.comps, b))


#= ADDITION =#

# addition of same-grade elements to produce homogeneous multivector
function +(as::HomogeneousMultivector{k}...) where k
	T = best_type(Multivector, as...; grade=Val(k))
	Σ = zero(T)
	for a ∈ as, u ∈ blades(a)
		add!(Σ, u)
	end
	Σ
	# sum(sum, (blades(u) for u in as))
end

# addition of possibly mixed grade elements to produce inhomogeneous multivector
function +(as::AbstractMultivector...)
	T = best_type(MixedMultivector, as...)
	Σ = zero(T)
	for a ∈ as, u ∈ blades(a)
		add!(Σ, u)
	end
	Σ
end

# addition of scalar to grade 0 objects
+(a::Blade{sig,0}, b::Scalar) where sig = Blade{sig,0}(a.coeff + b, a.ublade)
	
function +(a::Multivector{sig,0}, b::Scalar) where sig
	ab = zero(best_type(Multivector, a, grade=Val(0), el=typeof(b)))
	add!(ab, a)
	ab[] += b
	ab
end

function +(a::AbstractMultivector, b::Scalar)
	ab = zero(best_type(MixedMultivector, a, el=typeof(b)))
	add!(ab, a)
	ab[] += b
	ab
end

+(a::Scalar, b::AbstractMultivector) = +(b, a)


# additive inverses

-(a::Blade) = typeof(a)(-a.coeff, a.ublade)
-(a::AbstractMultivector) = -one(eltype(a))*a

-(a::AbstractMultivector, b::AbstractMultivector) = a + (-b)
-(a::AbstractMultivector, b::Scalar) = a + (-b)
-(a::Scalar, b::AbstractMultivector) = a + (-b)






#= ALGEBRAIC PRODUCTS =#

# TODO: Multivector{k}*Multivector{0} should give Multivector{k}

"""
	a*b
	geometric_prod(a, b)

Geometric product of multivectors `a` and `b`.

Returns a `MixedMultivector`, unless both arguments are `Blade`s, in which
case the product is a `Blade`.

Examples
===
```jldoctests
julia> @basis x y z
[ Info: Defined basis blades x, y, z, xy, xz, yz, xyz

julia> x*(1 + y)
MixedMultivector{⟨x+,y+,z+⟩, Array{Float64,1}}:
 1.0 x
 1.0 xy

julia> 2x*2x
Grade-0 Blade{⟨x+,y+,z+⟩, 0, Float64, UInt64}:
 4.0
```
"""
geometric_prod

function geometric_prod(a::Blade, b::Blade)
	sig = shared_sig(a, b)
	factor, u = ubladeprod(sig, a.ublade, b.ublade)
	Blade{sig}(factor*(a.coeff*b.coeff), u)
end
function geometric_prod(a::AbstractMultivector, b::AbstractMultivector)
	ab = zero(best_type(MixedMultivector, a, b))
	for u ∈ blades(a), v ∈ blades(b)
		add!(ab, geometric_prod(u, v))
	end
	ab
end

function homogeneous_prod(a::Blade, b::Blade, k)
	if ublade_grade(ublade_xor(a.ublade, b.ublade)) == k
		geometric_prod(a, b)
	else
		zero(best_type(Blade, a, b; grade=Val(k)))
	end
end

function homogeneous_prod(a::HomogeneousMultivector{k1}, b::HomogeneousMultivector{k2}, k) where {k1,k2}
	if k ∈ abs(k1 - k2):2:(k1 + k2)
		_homogeneous_prod(a, b, Val(k))
	else
		zero(best_type(Multivector, a, b; grade=Val(k)))
	end
end
function _homogeneous_prod(a::AbstractMultivector, b::AbstractMultivector, ::Val{k}) where k
	ab = zero(best_type(Multivector, a, b; grade=Val(k)))
	0 <= k <= dimension(a) || return ab
	for u ∈ blades(a), v ∈ blades(b)
		add!(ab, homogeneous_prod(u, v, k))
	end
	ab
end

function graded_prod(a::AbstractMultivector, b::AbstractMultivector, grade_selector)
	ab = zero(best_type(MixedMultivector, a, b))
	for u ∈ blades(a), v ∈ blades(b)
		k = grade_selector(grade(u), grade(v))
		0 <= k <= dimension(ab) || continue
		add!(ab, homogeneous_prod(u, v, k))
	end
	ab
end
graded_prod(a::HomogeneousMultivector, b::HomogeneousMultivector, grade_selector) =
	homogeneous_prod(a, b, grade_selector(grade(a), grade(b)))

"""
	*(::AbstractMultivector, ::AbstractMultivector)

Geometric product of multivectors. See also `GeometricAlgebra.geometric_prod`.
"""
*(a::AbstractMultivector, b::AbstractMultivector) = geometric_prod(a, b)


# faster versions for vector storage types
let T = CompositeMultivector{<:AbstractVector}
	@eval *(a::$T, b::$T) = geometric_prod_gen(a, b)
	@eval _homogeneous_prod(a::$T, b::$T, k) = homogeneous_prod_gen(a, b, k)
end



"""
```
a⋅b
dot(a, b)
```

"Fat" dot product of multivectors, giving the lowest-possible-grade part of the geometric product.
Defined by ``a⋅b = ∑_{ij} ⟨⟨a⟩_i * ⟨b⟩_j⟩_|i - j|``.
"""
dot(a, b) = graded_prod(a, b, abs∘-)
const ⋅ = dot

"""
```
a∧b
wedge(a, b)
```

Wedge product of multivectors, giving the highest-possible-grade part of the geometric product `a*b`.
Defined by ``a∧b = ∑_{ij} ⟨⟨a⟩_i * ⟨b⟩_j⟩_(i + j)``.
"""
wedge(a, b) = graded_prod(a, b, +)
const ∧ = wedge

## DEPRECIATED: just do `scalar(a*b)`
# """
# ```
# a∗b
# scalar_prod(a, b)
# ```
#
# Scalar part of multivector product, equivalent to `grade(a*b, 0)`.
# """
# scalar_prod(a, b) = homogeneous_prod(a, b, 0)
# const ∗ = scalar_prod

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




#= REVERSION, INVOLUTION, DUALITY, ETC =#

# @pure ?
reversionsign(k::Integer, c) = mod(k, 4) <= 1 ? c : -c
reversionsign(k::Integer) = reversionsign(k, 1)

reversion(a::Number) = a
reversion(a::HomogeneousMultivector{k}) where k = reversionsign(k, a)
reversion(a::MixedMultivector) = mapcomps(u -> reversion(u).coeff, a)
~(a::AbstractMultivector) = reversion(a)

involute(a::HomogeneousMultivector) = iseven(grade(a)) ? a : -a
involute(a::MixedMultivector) = mapcomps(u -> (iseven(grade(u)) ? u : -u).coeff, a)

# is this appropriate?
Base.adjoint(a::AbstractMultivector) = reversion(a)

# don't know what this operation is yet
conjugate(a::AbstractMultivector{sig}) where sig = mapcomps(u -> ublade_square(sig, u.ublade)*u.coeff, a)



"""
```
★(a)
hodgedual(a)
```

Hodge dual or 'star' operator, satisfying
``a∧★(b) == (a∗b)vol(a)``.(?)
"""
hodgedual(a::AbstractMultivector) = reversion(a)⋅vol(a)
# const ★ = hodgedual


dual(a::AbstractMultivector) = a⋅vol(a)


#= NORMS =#

"""
	rsqrt(x)

Real-valued square root defined on all real numbers as `sign(x)sqrt(abs(x))`.

Examples
===
```jldoctest; setup = :( using GeometricAlgebra: rsqrt )
julia> rsqrt(-9)
-3.0
```
"""
rsqrt(x) = sign(x)sqrt(abs(x))


Base.abs2(a::Blade) = abs2(a.coeff)*ublade_square(signature(a), a.ublade)
Base.abs2(a) = a'a

Base.abs(a) = let n = abs2(a)
	isscalar(n) || @warn "norm is not scalar" n
	rsqrt(scalar(n)) # is this appropriate?
end

norm_unit(a::Blade) = (a.coeff, oneunit(a))
norm_unit(a::AbstractMultivector) = let n = abs(a)
	(n, a/n)
end





#= MULTIPLICATIVE INVERSES =#

# blades have an inverse iff their square is nonzero
# homogeneous multivectors are invertible `sometimes` - TODO figure out when

# Base.inv(a::Blade) = a/scalar(a*a) # blades are guaranteed to have scalar squares
Base.inv(a::Blade) = a/(reversionsign(grade(a))ublade_square(signature(a), a.ublade)) # blades are guaranteed to have scalar squares
function Base.inv(a::AbstractMultivector)
	a² = a^2
	isscalar(a²) && return (a/scalar(a²))::typeof(a)

	inv_matrixmethod(a)
	# error("cannot find inverse of multivector $a")
end

# TODO use sparse matrices since they're 2^dimension × 2^dimension
function inv_matrixmethod(a::T) where T<:MixedMultivector{sig,C} where {sig,C<:AbstractVector}
	M = hcat([(a*Blade{sig}(one(eltype(a)), ublade)).comps for ublade ∈ 0:unsigned(2^dimension(a) - 1)]...)
	T(M\one(T).comps)
end

function inv_matrixmethod(a::CompositeMultivector{C}) where {C}
	sig = signature(a)
	d = 2^dimension(sig)
	M = Matrix{eltype(a)}(undef, d, d)
	for j ∈ 1:d
		u = Blade{sig}(one(eltype(a)), unsigned(j - 1))
		au = a*u
		for i ∈ 1:d
			M[i, j] = getcomp(au, unsigned(i - 1))
		end
	end
	x = zeros(eltype(a), d)
	x[1] = one(eltype(a))
	MixedMultivector{sig,C}(M\x)
end


/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
/(a::AbstractMultivector, b::Number) = a*inv(b)
/(a::Number, b::AbstractMultivector) = a*inv(b)

\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b
\(a::AbstractMultivector, b::Number) = inv(a)*b
\(a::Number, b::AbstractMultivector) = inv(a)*b





# EXPONENTIATION

function Base.literal_pow(::typeof(^), a::Blade{sig,k,T,B}, ::Val{2})::Blade{sig,0,T,B} where {sig,k,T,B}
	a.coeff^2*reversionsign(grade(a))ublade_square(signature(a), a.ublade)
end

function ^(a::Blade, p::Integer)
	p >= 0 || return inv(a)^abs(p)
	blade = iseven(p) ? one(a) : oneunit(a)
	square, _ = ubladeprod(signature(a), a.ublade, a.ublade)
	factor = square^fld(p, 2) # floored division
	coeff = factor*a.coeff^p
	coeff*blade
end

function powbysquaring(a::AbstractMultivector, p::Integer)
	p >= 0 || return powbysquaring(inv(a), abs(p))
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
^(a::AbstractMultivector, p) = error("unsupported multivector exponent type $(typeof(p))")


SERIES_ORDER = 10
set_order(n) = begin
	global SERIES_ORDER
	SERIES_ORDER = n
end
function eval_taylor_series(a::AbstractMultivector, func)
	T = eltype(a)
	order = SERIES_ORDER
	series = func(Taylor1(T, order))
	series(a)
end

function Base.exp(a::AbstractMultivector)
	a² = a^2
	if isscalar(a²)
		s = scalar(a²)
		if iszero(s)
			one(a) + a
		else
			coeff, â = norm_unit(a)
			s > 0 ? cosh(coeff) + sinh(coeff)â : cos(coeff) + sin(coeff)â
		end
	else
		# fallback to explicit computation
		eval_taylor_series(a, exp)
	end
end


# Use the identity ``log(x) = 2atanh((x - 1)/(x + 1))`` and hope for the best
Base.log(a::AbstractMultivector) = 2atanh((a - 1)/(a + 1))


# TRIGONOMETRIC FUNCTIONS

# even trig functions (e.g., cos) of blades produce scalars, since even powers are scalar
function apply_trig_even(a, sqsign, pos, neg)
	if iszero(sqsign)
		one(a)
	else
		coeff = abs(a)
		sqsign > 0 ? pos(coeff) : neg(coeff)
	end
end
# odd (e.g., sin) functions produce blades parallel to the original
function apply_trig_odd(a, sqsign, pos, neg)
	if iszero(sqsign)
		a
	else
		coeff, â = norm_unit(a)
		sqsign > 0 ? pos(coeff)*â : neg(coeff)*â
	end
end


for (trigfn,  (pos,    neg,    even )) ∈ [
	:cos   => (:cos,   :cosh,  true ),
	:cosh  => (:cosh,  :cos,   true ),
	:sin   => (:sin,   :sinh,  false),
	:sinh  => (:sinh,  :sin,   false),
	:tan   => (:tan,   :tanh,  false),
	:tanh  => (:tanh,  :tan,   false),
	:asin  => (:asin,  :asinh, false),
	:asinh => (:asinh, :asin,  false),
	:atan  => (:atan,  :atanh, false),
	:atanh => (:atanh, :atan,  false),
]
	apply_trig = even ? apply_trig_even : apply_trig_odd
	@eval function Base.$trigfn(a::AbstractMultivector)
		a² = a*a
		if isscalar(a²)
			$apply_trig(a, scalar(a²), $pos, $neg)
		else
			eval_taylor_series(a, $trigfn)
		end
	end
end

for (trigfn, invfn) ∈ [
	:sec => :cos,
	:sech => :cosh,
	:csc => :sin,
	:csch => :sinh,
	:cot => :tan,
	:coth => :tanh,
]
	@eval Base.$trigfn(a::AbstractMultivector) = inv(Base.$invfn(a))
end