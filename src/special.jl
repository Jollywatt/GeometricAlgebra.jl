function eval_evenodd_trig(a, a², pos, neg, ::Val{even}) where even
	s = sign(a²)
	norm = sqrt(abs(a²))
	if even
		s > 0 ? pos(norm) : s < 0 ? neg(norm) : one(a)
	else
		s > 0 ? pos(norm)*a/norm : s < 0 ? neg(norm)*a/norm : a
	end
end

function eval_series(a, fn, order=30)
	T = eltype(a)
	series = fn(Taylor1(T, order))
	series(a)
end

function Base.exp(a::AbstractMultivector)
	a² = a*a
	if isscalar(a²)
		s = scalar(a²)
		cosha = eval_evenodd_trig(a, s, cosh, cos, Val(true))
		sinha = eval_evenodd_trig(a, s, sinh, sin, Val(false))
		cosha + sinha
	else
		eval_series(a, exp)
	end
end


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