@static if VERSION < v"1.7"
	ismutabletype(a::DataType) = a.mutable
end

#= Array-like Type Utilities =#


# length of fixed-size contained from its type
length_from_type(::Type{<:NTuple{N}}) where {N} = (N)
length_from_type(T::Type{<:StaticVector}) = (length(T))

with_eltype(::Type{<:Vector}, T) = Vector{T}
with_eltype(::Type{<:MVector{N}}, T) where {N} = MVector{N,T}
with_eltype(::Type{<:SVector{N}}, T) where {N} = SVector{N,T}
with_eltype(::Type{<:NTuple{N}}, T) where {N} = NTuple{N,T}



zeroslike(::Type{Vector{T}}, n) where {T<:Number} = zeros(T, n)
zeroslike(::Type{<:MVector{N,T}}, n) where {N,T} = zeros(MVector{n,T})
zeroslike(::Type{<:SVector{N,T}}, n) where {N,T} = zeros(SVector{n,T})

oneslike(::Type{Vector{T}}, n) where {T<:Number} = ones(T, n)
oneslike(::Type{<:MVector{N,T}}, n) where {N,T} = ones(MVector{n,T})
oneslike(::Type{<:SVector{N,T}}, n) where {N,T} = ones(SVector{n,T})


#= Tools to handle symbolic eltypes =#

# Multivectors with purely symbolic components have eltype T<:Symbolic,
# but as soon as they are mixed with e.g., Int, have eltype Any.
# So we do want so support Vector{Any} storage types for CompositeMultivectors,
# but this should only happen for symbolic stuff

isrealzero(a::Number) = iszero(a)
isrealzero(a) = false
isrealone(a::Number) = isone(a)
isrealone(a) = false

realzero(T::Type{<:Number}) = zero(T)
realzero(::Type) = zero(Int)
realone(T::Type{<:Number}) = one(T)
realone(::Type) = one(Int)

oneslike(::Type{<:Vector}, dims...) = let a = ones(Int, dims...)
	convert(with_eltype(typeof(a), Any), a)
end
zeroslike(::Type{<:Vector}, dims...) = let a = zeros(Int, dims...)
	convert(with_eltype(typeof(a), Any), a)
end

numberorany(T::Type{<:Number}) = T
numberorany(::Type) = Any


#= Miscellaneous =#

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



function copy_setindex(a, val, I...)
	T = promote_type(eltype(a), typeof(val))
	a′ = with_eltype(typeof(a), T)(a)
	if ismutable(a)
		setindex!(a′, val, I...)
		a′
	else
		setindex(a′, val, I...)
	end
end