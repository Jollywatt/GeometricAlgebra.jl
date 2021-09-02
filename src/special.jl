function eval_evenodd_trig(a, a², pos, neg, ::Val{even}) where even
	s = sign(a²)
	norm = sqrt(abs(a²))
	if even
		s > 0 ? pos(norm) : s < 0 ? neg(norm) : one(a)
	else
		s > 0 ? pos(norm)*a/norm : s < 0 ? neg(norm)*a/norm : a
	end
end

SERIES_ORDER = [150]

function eval_series(a, fn, order=SERIES_ORDER[])
	T = eltype(a)
	series = fn(Taylor1(T, order))
	series(a)
end

isapproxzero(x) = abs(x) <= eps(float(one(x)))/2

crudenorm(a::Scalar) = a
crudenorm(a::Blade) = a.coeff
crudenorm(a::CompositeMultivector) = sum(abs.(a.components))



function exp_with_scalar_square(a, a²)
	isapproxzero(a²) && return one(a) + a

	norm = sqrt(abs(a²))
	if a² > 0
		cosh(norm) + sinh(norm)/norm*a
	else
		cos(norm) + sin(norm)/norm*a
	end
end

function exp_series(a)
	# Series convergence is better when `a` is not too large.
	# Use the fact that ``exp(a) = exp(a/p)^p`` and choose `p`
	# so that the norm of `a` is of order one.

	norm = crudenorm(a) # always positive
	p = 2^max(0, round(Int, log2(norm))) # power of 2 for fast exponentiation
	a /= p

	result = one(a)
	term = one(a)

	# empirically chosen so that `term` usually underflows before `max_iters` is reached
	max_iters = ceil(Int, precision(eltype(a))^(3//4))
	for i = 1:max_iters
		term *= a/i
		if isapproxzero(crudenorm(term))
			break
		end
		result += term
	end

	result^p
end

function Base.exp(a::AbstractMultivector)
	# use fact that exp(scalar + a) = exp(scalar)exp(a)
	s = scalar(a)
	if !iszero(s)
		prefactor = exp(s)
		a -= s
	else
		prefactor = one(s)
	end

	a² = a*a
	if isscalar(a²)
		prefactor*exp_with_scalar_square(a, scalar(a²))
	else
		prefactor*exp_series(a - s)
	end
end

Base.exp(a::Blade) = exp_with_scalar_square(a, scalar(a*a))





for (fn,     pos,    neg,    even) ∈ [
	(:cosh,  :cosh,  :cos,   true),
	(:cos,   :cos,   :cosh,  true),
	(:sinh,  :sinh,  :sin,   false),
	(:sin,   :sin,   :sinh,  false),
	(:tanh,  :tanh,  :tan,   false),
	(:tan,   :tan,   :tanh,  false),
	(:asinh, :asinh, :asin,  false),
	(:asin,  :asin,  :asinh, false),
	(:atanh, :atanh, :atan,  false),
	(:atan,  :atan,  :atanh, false),
]
	@eval function Base.$fn(a::AbstractMultivector)
		a² = a*a
		if isscalar(a²)
			eval_evenodd_trig(a, scalar(a²), $pos, $neg, Val($even))
		else
			eval_series(a, $fn)
		end
	end
end

for (fn,    invfn) ∈ [
	(:sec,  :cos),
	(:sech, :cosh),
	(:csc,  :sin),
	(:csch, :sinh),
	(:cot,  :tan),
	(:coth, :tanh),
]
	@eval Base.$fn(a::AbstractMultivector) = inv(Base.$invfn(a))
end


# log(x) = 2atanh((x - 1)/(x + 1))
Base.log(a::AbstractMultivector) = 2atanh((a - 1)/(a + 1))
