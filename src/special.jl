#= Multiplicative Inverses =#

"""
	vector_repr(a::AbstractMultivector)

Vector representation of a multivector.

This is a linear operator. See also [`matrix_repr`](@ref).

# Examples
```jldoctest
julia> @basis 2
[ Info: Defined basis blades v, v1, v2, v12

julia> vector_repr(v1*v2) == matrix_repr(v1)vector_repr(v2)
true
```
"""
vector_repr(a::AbstractMultivector) = collect(MixedMultivector(a).comps)

"""
	matrix_repr(a::AbstractMultivector)

Matrix representation of a multivector.

This is an injective homomorphism from the geometric algebra
to a matrix subalgebra (i.e., it is linear, and preserves algebraic products).

See also [`vector_repr`](@ref).

# Examples
```jldoctest
julia> @basis 2
[ Info: Defined basis blades v, v1, v2, v12

julia> matrix_repr(v1)
4×4 Matrix{Int64}:
 0  1  0  0
 1  0  0  0
 0  0  0  1
 0  0  1  0

julia> matrix_repr(v1*v2) == matrix_repr(v1)matrix_repr(v2)
true
```
"""
matrix_repr(a::HomogeneousMultivector) = matrix_repr(MixedMultivector(a))
function matrix_repr(a::MixedMultivector)
	N, T = ncomponents(a), eltype(a)
	mat = Matrix{numberorany(T)}(undef, N, N)
	fill!(mat, numberzero(T))
	for (i, b) ∈ enumerate(basis(signature(a), grade=:all))
		mat[:,i] = MixedMultivector(a*b).comps
	end
	mat
end

function via_matrix_repr(f::Function, a::AbstractMultivector)
	m = matrix_repr(a)
	m′ = f(m)
	MixedMultivector{signature(a)}(m′[:,1])
end



function inv_matrix_method(a::CompositeMultivector)
	A = matrix_repr(a)
	id = MixedMultivector(one(a)).comps
	A⁻¹ = A\id
	MixedMultivector{signature(a)}(A⁻¹)
end

function inv_formula_method(a::CompositeMultivector)
	# In low dimensions, explicit formulae exist.
	# See https://doi.org/10.1016/j.amc.2017.05.027

	dim = dimension(a)
	ā = clifford_conj(a)

	if dim <= 2
		return ā/scalar(a*ā)

	elseif dim == 3
		āâã = ā*involution(a)*reversion(a)
		return āâã/scalar(a*āâã)

	elseif dim == 4
		b = ā*graded_multiply(a*ā) do k
			k ∈ (3, 4) ? -1 : 1
		end
		return b/scalar(a*b)

	elseif dim == 5
		b = ā*reversion(a*ā)
		c = b*graded_multiply(a*b) do k
			k ∈ (1, 4) ? -1 : 1
		end
		return c/scalar(a*c)

	else	
		throw("only implemented for dimensions 0:5")
	end
end

Base.inv(a::Blade) = a/scalar(a^2)
function Base.inv(a::CompositeMultivector)
	if dimension(a) <= 5
		inv_formula_method(a)
	else
		inv_matrix_method(a) # ?? replace with via_matrix_repr(inv, a)
	end
end

Base.:/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b

Base.:/(a::Scalar, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::Scalar) = inv(a)*b



#= Exponential =#

function exp_with_scalar_square(a, a²::Scalar)
	norm = sqrt(abs(a²))
	if iszero(norm)
		one(a) + a
	elseif a² > 0
		cosh(norm) + sinh(norm)/norm*a
	else
		cos(norm) + sin(norm)/norm*a
	end
end


infnorm(a::Blade) = abs(a.coeff)
infnorm(a::CompositeMultivector) = maximum(abs.(a.comps))

twonorm(a::Blade) = abs(a.coeff)
twonorm(a::CompositeMultivector) = sqrt(sum(abs2.(a.comps)))

exp_series(a::HomogeneousMultivector) = exp_series(MixedMultivector(a))
function exp_series(a::MixedMultivector{Sig,C}) where {Sig,C}
	# Series convergence is better when `a` is not too large.
	# Use the fact that ``exp(λa) = exp(λ)exp(a/p)^p`` and choose `p`
	# so that the 2-norm of `a` is of order one.

	λ = scalar(a)
	a -= λ

	norm = twonorm(a)
	p = 2^max(0, floor(Int, log2(norm))) # power of 2 for fast exponentiation
	a /= p

	term = one(a)
	result = copy(term)

	max_iters = 200
	for i in 1:max_iters
		term *= a/i
		infnorm(term) < eps(eltype(term)) && break
		result += term
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

Base.sqrt(a::Blade{Sig,0}) where {Sig} = Blade{Sig,0}(0 => sqrt(a.coeff))
Base.sqrt(a::CompositeMultivector) = isscalar(a) ? sqrt(scalar(a))one(a) : via_matrix_repr(sqrt, a)

Base.log(a::AbstractMultivector) = via_matrix_repr(log, a)



#= Trigonometric =#

function eval_evenodd_trig(a, a², pos, neg, ::Val{even}) where even
	s = sign(a²)
	norm = sqrt(abs(a²))
	if even
		s > 0 ? pos(norm) : s < 0 ? neg(norm) : one(a)
	else
		s > 0 ? pos(norm)*a/norm : s < 0 ? neg(norm)*a/norm : a
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

