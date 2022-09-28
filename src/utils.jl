@static if VERSION < v"1.7"
	ismutabletype(a::DataType) = a.mutable
end

#= Array-like Type Utilities =#

issetindexable(T::Type) = ismutabletype(T)
issetindexable(::Type{<:AbstractSparseArray}) = true
issetindexable(a) = issetindexable(typeof(a))
issetindexable(::Type{<:MVector{N,T}}) where {N,T} = isbitstype(T) # see https://github.com/JuliaArrays/StaticArrays.jl/issues/27

zeroslike(::Type{<:Array{T}}, dims...) where {T<:Number} = zeros(T, dims...)
zeroslike(::Type{<:MArray{N,T}}, dims...) where {N,T} = zeros(MArray{Tuple{dims...},T})
zeroslike(::Type{<:SArray{N,T}}, dims...) where {N,T} = zeros(SArray{Tuple{dims...},T})
zeroslike(::Type{<:SparseVector{T}}, dims...) where {T} = spzeros(T, dims...)

oneslike(::Type{<:Array{T}}, dims...) where {T<:Number} = ones(T, dims...)
oneslike(::Type{<:MArray{N,T}}, dims...) where {N,T} = ones(MArray{Tuple{dims...},T})
oneslike(::Type{<:SArray{N,T}}, dims...) where {N,T} = ones(SArray{Tuple{dims...},T})
oneslike(::Type{<:SparseVector{N,T}}, dims...) where {N,T} = sparse(ones(T, dims...))


with_eltype(::Type{<:Vector}, T) = Vector{T}
with_eltype(::Type{<:MVector{N}}, T) where {N} = MVector{N,T}
with_eltype(::Type{<:SVector{N}}, T) where {N} = SVector{N,T}
with_eltype(::Type{<:SparseVector}, T) = SparseVector{T}

function copy_setindex(a, val, I...)
	T = promote_type(eltype(a), typeof(val))
	a′ = with_eltype(typeof(a), T)(a)
	if issetindexable(a)
		setindex!(a′, val, I...)
		a′
	else
		setindex(a′, val, I...)
	end
end


#= Tools to handle symbolic eltypes =#

# Multivectors with purely symbolic components have eltype T<:Symbolic,
# but as soon as they are mixed with e.g., Int, have eltype Any.
# So we do want so support Vector{Any} storage types for CompositeMultivectors,
# but this should only happen for symbolic stuff

isnumberzero(a::Number) = iszero(a)
isnumberzero(a) = false
isnumberone(a::Number) = isone(a)
isnumberone(a) = false

numberzero(T::Type{<:Number}) = zero(T)
numberzero(::Type) = zero(Int)
numberone(T::Type{<:Number}) = one(T)
numberone(::Type) = one(Int)

oneslike(T::Type{<:Array}, dims...) = convert(Array{Any}, ones(Int, dims...))
zeroslike(T::Type{<:Array}, dims...) = convert(Array{Any}, zeros(Int, dims...))

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



