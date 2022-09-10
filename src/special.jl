#= Multiplicative Inverses =#

Base.inv(a::Blade) = a/scalarpart(a^2)
function Base.inv(a::CompositeMultivector)
	aã = a*reversion(a)
	isscalar(aã) && return reversion(a)/scalarpart(aã)

	a² = a*a
	isscalar(a²) && return a/scalarpart(a²)
	
	throw("not yet implemented for multivectors for which `A*reversion(A)` is non-scalar")
end

Base.:/(a::AbstractMultivector, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::AbstractMultivector) = inv(a)*b

Base.:/(a::Number, b::AbstractMultivector) = a*inv(b)
Base.:\(a::AbstractMultivector, b::Number) = inv(a)*b

function exp_with_scalar_square(a, a²::Number)
	norm = sqrt(abs(a²))
	if a² > 0
		cosh(norm) + sinh(norm)/norm*a
	else
		cos(norm) + sin(norm)/norm*a
	end
end


function Base.exp(a::AbstractMultivector)
	a² = a^2
	if isscalar(a²)
		exp_with_scalar_square(a, scalarpart(a²))
	else
		throw("not yet implemented for multivectors with non-scalar square")
	end
end