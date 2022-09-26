#= Multiplicative Inverses =#

matrix_repr(a::HomogeneousMultivector) = matrix_repr(MixedMultivector(a))
function matrix_repr(a::MixedMultivector)
	N = ncomponents(a)
	if eltype(a) <: Number
		mat = zeros(eltype(a), N, N)
	else
		mat = Matrix{Any}(undef, N, N)
		fill!(mat, 0)
	end
	for (i, b) ∈ enumerate(basis(signature(a), grade=:all))
		mat[:,i] = MixedMultivector(a*b).comps
	end
	mat
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
		return ā/scalarpart(a*ā)

	elseif dim == 3
		āâã = ā*involution(a)*reversion(a)
		return āâã/scalarpart(a*āâã)

	elseif dim == 4
		b = ā*graded_multiply(a*ā) do k
			k ∈ (3, 4) ? -1 : 1
		end
		return b/scalarpart(a*b)

	elseif dim == 5
		b = ā*reversion(a*ā)
		c = b*graded_multiply(a*b) do k
			k ∈ (1, 4) ? -1 : 1
		end
		return c/scalarpart(a*c)

	else	
		throw("only implemented for dimensions 0:5")
	end
end

Base.inv(a::Blade) = a/scalarpart(a^2)
function Base.inv(a::CompositeMultivector)
	if dimension(a) <= 5
		inv_formula_method(a)
	else
		inv_matrix_method(a)
	end
end

Base.:/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b

Base.:/(a::Scalar, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::Scalar) = inv(a)*b


#= Exponential =#

function exp_with_scalar_square(a, a²::Scalar)
	norm = sqrt(abs(a²))
	if a² > 0
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
	# Use the fact that ``exp(a) = exp(a/p)^p`` and choose `p`
	# so that the 2-norm of `a` is of order one.
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

	result^p
end


function Base.exp(a::AbstractMultivector)
	a² = a^2
	if isscalar(a²)
		exp_with_scalar_square(a, scalarpart(a²))
	else
		exp_series(a)
	end
end
