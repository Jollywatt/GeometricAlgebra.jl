#= Multiplicative Inverses =#

"""
	matrix_repr(a::AbstractMultivector, k=0:dim)

Matrix representation of the grade `k` parts of a multivector.

By default, the full ``2^d × 2^d`` linear representation is used in ``d`` dimensions.
Smaller representations can be used for elements in

- the even subalgebra, `k=0:2:dim`
- the scalar-pseudoscalar subalgebra, `k=(0, dim)`

by restricting `k` to those grades.

# Examples
```jldoctest
julia> @basis 2
[ Info: Defined basis blades v1, v2, v12, I in Main

julia> matrix_repr(1 + 7v12)
4×4 Matrix{Int64}:
 1   0  0  -7
 0   1  7   0
 0  -7  1   0
 7   0  0   1

julia> matrix_repr(1 + 7v12, (0, 2))
2×2 Matrix{Int64}:
 1  -7
 7   1

julia> matrix_repr(v1*v2) == matrix_repr(v1)matrix_repr(v2)
true
```
"""
function matrix_repr(a::AbstractMultivector, k=0:dimension(a))
	a = grade(a, k)
	N, T = ncomponents(a), eltype(a)
	mat = zeroslike(Matrix{eltype(a)}, N, N)
	M = Multivector{signature(a),grade(a)}
	for (i, b) ∈ enumerate(basis(M))
		mat[:,i] = (a*b).comps
	end
	mat
end

@symbolic_optim function matrix_repr(a::Multivector, ::Val{K}) where {K}
	matrix_repr(a, K)
end

function via_matrix_repr(f::Function, a::AbstractMultivector)
	k = resulting_grades(Val(:subalgebra), dimension(a), grade(a))
	m = matrix_repr(a, Val(k))
	m′ = f(m)
	T = componentstype(signature(a), size(m′, 1), eltype(m′))
	Multivector{signature(a),k}(convert(T, m′[:,1]))
end


function inv_matrix_method(a::Multivector)
	k = resulting_grades(Val(:subalgebra), dimension(a), grade(a))
	A = matrix_repr(a, k)
	id = zeros(size(A, 1))
	id[1] = 1
	A⁻¹ = A\id
	Multivector{signature(a),k}(A⁻¹)
end

function inv_formula_method(a::Multivector)
	# In low dimensions, explicit formulae exist.
	# See https://doi.org/10.1016/j.amc.2017.05.027

	dim = dimension(a)
	scalar_square = grade(a) ∈ (0, 1, dim - 1, dim) # (pseudo)(scalars|vectors) always square to scalars
	scalar_square && return a/scalar(a^2)

	ā = clifford_conj(a)
	aā = a*ā

	if dim <= 2 || isscalar(aā)
		return ā/scalar(aā)

	elseif dim == 3
		āâã = ā*involution(a)*reversion(a)
		return āâã/scalar(a*āâã)

	elseif dim == 4
		b = ā*graded_multiply(aā) do k
			k ∈ (3, 4) ? -1 : 1
		end
		return b/scalar(a*b)

	elseif dim == 5
		b = ā*reversion(aā)
		c = b*graded_multiply(a*b) do k
			k ∈ (1, 4) ? -1 : 1
		end
		return c/scalar(a*c)

	else	
		inv_matrix_method(a) # ?? replace with via_matrix_repr(inv, a)
	end
end

Base.inv(a::BasisBlade) = a/scalar(a^2)
Base.inv(a::Multivector) = inv_formula_method(a)

Base.:/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b

Base.:/(a::Scalar, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::Scalar) = inv(a)*b



#= Exponential =#

function exp_with_scalar_square(a, a²::Scalar)
	iszero(a²) && return one(a)/one(a²) + a
	if a² isa Real && a² < 0
		norm = sqrt(-a²)
		cos(norm) + sin(norm)/norm*a
	else
		norm = sqrt(a²)
		cosh(norm) + sinh(norm)/norm*a
	end
end


infnorm(a::BasisBlade) = abs(a.coeff)
infnorm(a::Multivector) = maximum(abs.(a.comps))

twonorm(a::BasisBlade) = abs(a.coeff)
twonorm(a::Multivector) = sqrt(sum(abs2.(a.comps)))


function exp_series(a::Multivector)
	# Series convergence is better when `a` is not too large.
	# Use the fact that `exp(a) = exp(a/p)^p` and choose `p`
	# so that the 2-norm of `a` is of order one.

	# Also use the fact that `exp(λ + a) = exp(λ)exp(a)` for scalar `λ`
	λ = scalar(a)
	a -= λ

	norm = twonorm(a)
	p = 2^floor(Int, log2(max(1, norm))) # power of 2 for fast exponentiation
	a /= p

	term = one(a)
	result = grade(term, resulting_grades(Val(:subalgebra), dimension(a), grade(a)))

	max_iters = 200
	for i in 1:max_iters
		term *= a/i
		infnorm(term) < eps(real(eltype(term))) && break
		result = add!(result, term)
	end

	exp(λ)result^p
end


function Base.exp(a::AbstractMultivector)
	a² = a^2
	if isscalar(a²)
		exp_with_scalar_square(a, scalar(a²))
	else
		exp_series(a)
	end
end



#= Roots and Logs =#

# TODO: sqrt_formula_method?
function Base.sqrt(a::AbstractMultivector)

	if isscalar(a)
		s = scalar(a)
		s >= 0 && return sqrt(s)one(a)
		if pseudoscalar_square(a) < 0
			sqrt(abs(s))unit_pseudoscalar(a)
		end
	end

	a² = a^2
	if isscalar(a²)
		s = scalar(a²)
		λ = sqrt(abs(s))
		if s < 0
			return (a + λ)/sqrt(2λ)
		elseif s > 0
			if pseudoscalar_square(a) < 0
				I = unit_pseudoscalar(a)
				return (a + λ*I)/(1 + I)sqrt(λ)
			end
		end
	end

	via_matrix_repr(sqrt, a)
end

Base.log(a::AbstractMultivector) = via_matrix_repr(log, a)



#= Trigonometric =#

function eval_evenodd_trig(a, a², pos, neg, ::Val{even}) where even
	iszero(a²) && return even ? one(a) : a
	if a² isa Real && a² < 0
		norm = sqrt(-a²)
		even ? neg(norm) : neg(norm)*a/norm
	else
		norm = sqrt(a²)
		even ? pos(norm) : pos(norm)*a/norm
	end
end


for (pos,    neg,    even) ∈ [
	(:cosh,  :cos,   true),
	(:cos,   :cosh,  true),
	(:sinh,  :sin,   false),
	(:sin,   :sinh,  false),
	(:tanh,  :tan,   false),
	(:tan,   :tanh,  false),
	(:asinh, :asin,  false),
	(:asin,  :asinh, false),
	(:atanh, :atan,  false),
	(:atan,  :atanh, false),
]
	@eval function Base.$pos(a::AbstractMultivector)
		a² = a*a
		if isscalar(a²)
			eval_evenodd_trig(a, scalar(a²), $pos, $neg, Val($even))
		else
			via_matrix_repr($pos, a)
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

