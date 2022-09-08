promote_to(T, x) = convert(promote_type(T, typeof(x)), x)

zeroslike(::Type{Vector{T}}, n) where {T} = zeros(T, n)

with_eltype(::Type{<:Vector}, T) = Vector{T}
# with_size(V::Type{<:Vector}, N) = V

shared_sig(::OrType{<:AbstractMultivector{Sig}}...) where {Sig} = true
shared_sig(::OrType{<:AbstractMultivector}...) = false

function __init__()
	if isdefined(Base.Experimental, :register_error_hint)
		Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
			# try to detect cases where the method error is due to mixing different metric signatures
			if all(isa.(argtypes, Type{<:AbstractMultivector})) && !shared_sig(argtypes...)
				println(io, "\nPerhaps the multivectors have incompatible metric signatures?")
			end
		end
	end
end