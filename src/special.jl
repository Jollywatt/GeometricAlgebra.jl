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

"""
	try_ensure_real(a::Multivector)

Tries to convert a complex-valued `Multivector` into an equivalent real one,
returning the original multivector if it failed.

If an algebra ``G`` has a commuting pseudoscalar squaring to ``-1``, then there
is a canonical map ``G ⊗ ℂ ↦ G`` from the complexified algebra into itself given
my sending the imaginary unit to the pseudoscalar.
"""
function try_ensure_real(a::Multivector)
	eltype(a) <: Real && return a
	pseudoscalar_square(a) < 0 && isodd(dimension(a)) || return a
	M = constructor(a)
	M(real.(a.comps)) + M(imag.(a.comps))*unit_pseudoscalar(a)
end

function via_matrix_repr(f::Function, a::AbstractMultivector)
	wasreal = eltype(a) <: Real
	k = resulting_grades(Val(:subalgebra), dimension(a), grade(a))
	m = matrix_repr(a, Val(k))
	m′ = f(m)
	T = componentstype(signature(a), size(m′, 1), eltype(m′))
	a′ = Multivector{signature(a),k}(convert(T, m′[:,1]))
	wasreal ? try_ensure_real(a′) : a′
end


function inv_matrix_method(a::Multivector)
	k = resulting_grades(Val(:subalgebra), dimension(a), grade(a))
	A = matrix_repr(a, k)
	id = zeros(size(A, 1))
	id[1] = 1
	A⁻¹ = A\id
	Multivector{signature(a),k}(A⁻¹)
end

function has_scalar_square(::T) where T<:AbstractMultivector
	dim = dimension(T)
	grade(T) ∈ (0, 1, dim - 1, dim)
end

function inv_formula_method(a::Multivector{Sig,K}) where {Sig,K}
	# In low dimensions, explicit formulae exist.
	# See https://doi.org/10.1016/j.amc.2017.05.027

	has_scalar_square(a) && return a/scalar(a^2)

	dim = dimension(a)

	ā = clifford_conj(a)
	aā = a*ā

	if dim <= 2 || isscalar(aā)
		return ā/scalar(aā)

	elseif dim == 3
		āâã = graded_prod(Val(K), ā, reversion(aā))
		return āâã/scalar_prod(a, āâã)

	elseif dim == 4
		b = graded_multiply(aā) do k
			k ∈ (3, 4) ? -1 : 1
		end
		c = graded_prod(Val(K), ā, b)
		return c/scalar_prod(a, c)

	elseif dim == 5
		b = ā*reversion(aā)
		c = graded_multiply(a*b) do k
			k ∈ (1, 4) ? -1 : 1
		end
		d = graded_prod(Val(K), b, c)
		return d/scalar_prod(a, d)

	else
		mask = nonzerobitmask(a)
		effective_dim = count_ones(mask)
		if effective_dim <= 5
			via_subalgebra_mask(inv, a, mask)
		elseif effective_dim < dim
			via_subalgebra_mask(inv_matrix_method, a, mask)
		else
			inv_matrix_method(a)
		end
	end
end

function max_nonzero_grade(a::Multivector) 
	i = findlast(!iszero, a.comps)
	isnothing(i) && return 0
	count_ones(componentbits(a)[i])
end

"""
	inv_flv_method(::AbstractMultivector)

Inverse of a multivector using the Faddeev–LeVerrier algorithm [^1].

This algorithm requires ``2^{d - 1}`` many geometric multiplications, where ``d`` is
the dimension of the algebra.

[^1]: "Algorithmic Computation of Multivector Inverses and Characteristic Polynomials in Non-degenerate Clifford Algebras", [Dimiter2024](@cite).
"""
function inv_flv_method(a::Multivector)
	n = 2^(cld(max_nonzero_grade(a), 2) + 1)
	m = one(a)
	am = a
	c = -n*scalar(a)
	for k = 2:n
		m = am + c
		am = a*m
		c = -n/k*scalar(am)
	end
	iszero(c) && error("Multivector has no inverse")
	-m/c
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
			if pseudoscalar_square(a) < 0 && iseven(a)
				I = unit_pseudoscalar(a)
				return (a + λ*I)/(1 + I)sqrt(λ)
			end
		end
	end


	# try dropping any unused basis vectors
	# even reducing the total dimension by one seems to pay off
	mask = nonzerobitmask(a)
	if count_ones(mask) < dimension(a)
		return via_subalgebra_mask(sqrt, a, mask)
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

function project_to_subalgebra(a::Multivector, mask::Unsigned)
	allbits = componentbits(a)

	indices = bits_to_indices(mask)
	subsig = Tuple(basis_vector_square.(Ref(signature(a)), indices))

	if all(==(1), subsig)
		subsig = length(subsig)
	end

	subsig

	k = promote_grades(dimension(subsig), grade(a))
	a′ = Multivector{subsig,k}(a.comps[allbits .& mask .== allbits])

end

function embed_in_superalgebra(sig, a::Multivector, mask::Unsigned)
	k = promote_grades(dimension(sig), grade(a)...)
	T = Multivector{sig,k}
	allbits = componentbits(T)

	j = 1
	b = zeroslike(a.comps, length(allbits))
	for (i, bits) ∈ enumerate(allbits)
		if bits & mask == bits
			b[i] = a.comps[j]
			j += 1
		end
	end
	T(b)
end

function nonzerobitmask(a::Multivector)
	usedbits = zero(UInt)
	for (_, bits) in nonzero_components(a)
		usedbits |= bits
	end
	usedbits
end


function via_subalgebra_mask(fn, a::Multivector, mask)
	a′ = project_to_subalgebra(a, mask)
	b = fn(a′)
	embed_in_superalgebra(signature(a), b, mask)
end