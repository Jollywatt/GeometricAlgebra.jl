promote_to(T, x) = convert(promote_type(T, typeof(x)), x)


zeroslike(::Type{Vector{T}}, n) where {T} = zeros(T, n)
zeroslike(::Type{<:MVector{N,T}}, n) where {N,T} = zeros(MVector{n,T})

oneslike(::Type{Vector{T}}, n) where {T} = ones(T, n)
onesslike(::Type{<:MVector{N,T}}, n) where {N,T} = ones(MVector{n,T})

with_eltype(::Type{<:Vector}, T) = Vector{T}
with_eltype(::Type{<:MVector{N}}, T) where N = MVector{N,T}
# with_size(V::Type{<:Vector}, N) = V

shared_sig(::OrType{<:AbstractMultivector{Sig}}...) where {Sig} = true
shared_sig(::OrType{<:AbstractMultivector}...) = false

function __init__()
	if isdefined(Base.Experimental, :register_error_hint)
		Base.Experimental.register_error_hint(MethodError) do io, err, argtypes, kwargs
			# try to detect cases where the method error is due to mixing different metric signatures
			if all(isa.(argtypes, Type{<:AbstractMultivector})) && !shared_sig(argtypes...)
				println(io, "\nPerhaps the multivectors have incompatible metric signatures?")
			end 
		end
	end
end

const SUBSCRIPT_DIGITS   = collect("₀₁₂₃₄₅₆₇₈₉")
const SUPERSCRIPT_DIGITS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(n::Integer)   = join(SUBSCRIPT_DIGITS[begin .+ reverse(digits(n))])
superscript(n::Integer) = join(SUPERSCRIPT_DIGITS[begin .+ reverse(digits(n))])