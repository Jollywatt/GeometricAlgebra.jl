#= Tools to handle symbolic eltypes =#

# Multivectors with purely symbolic components have SymbolicUtils.Symbolic eltype.
# For sym::SymbolicUtils.Sym, iszero(sym) is a symbolic expression “sym == 0”,
# but we generally want `iszero(sym) == false` because it *could* be non-zero.
# Furthermore, If a symbolic Blade is promoted to a KVector, its eltype will
# be Union{Int,Sym} or Any, so we want to handle Any eltypes gracefully.


isnumberzero(a::Number) = iszero(a)
isnumberone(a::Number) = isone(a)

# don't treat a symbolic type as one/zero, even if it could represent one
isnumberzero(a) = false
isnumberone(a) = false

numberzero(T::Type{<:Number}) = zero(T)
numberone(T::Type{<:Number}) = one(T)

numberzero(::Type) = zero(Int)
numberone(::Type) = one(Int)

numberorany(T::Type{<:Number}) = T
numberorany(::Type) = Any


#= Array-like Type Utilities =#

# initialize zeros/ones while abstracting away array type
zeroslike(::Type{<:Array{T}}, dims...) where {T<:Number} = zeros(T, dims...)
 oneslike(::Type{<:Array{T}}, dims...) where {T<:Number} =  ones(T, dims...)
zeroslike(::Type{<:Array}, dims...) = convert(Array{Any}, zeros(Int, dims...))
 oneslike(::Type{<:Array}, dims...) = convert(Array{Any},  ones(Int, dims...))
zeroslike(::Type{<:MArray{N,T}}, dims...) where {N,T} = zeros(MArray{Tuple{dims...},numberorany(T)})
 oneslike(::Type{<:MArray{N,T}}, dims...) where {N,T} =  ones(MArray{Tuple{dims...},numberorany(T)})
zeroslike(::Type{<:SArray{N,T}}, dims...) where {N,T} = zeros(SArray{Tuple{dims...},numberorany(T)})
 oneslike(::Type{<:SArray{N,T}}, dims...) where {N,T} =  ones(SArray{Tuple{dims...},numberorany(T)})
zeroslike(::Type{<:SparseVector{Tv}}, dims...) where {Tv} = spzeros(Tv, dims...)
 oneslike(::Type{<:SparseVector{Tv}}, dims...) where {Tv} = sparse(ones(Tv, dims...))



@static if VERSION < v"1.7"
	ismutabletype(a::DataType) = a.mutable
end

# whether array supports `setindex!`
issetindexable(T::Type) = ismutabletype(T)
issetindexable(::Type{<:AbstractSparseArray}) = true
issetindexable(a) = issetindexable(typeof(a))
issetindexable(::Type{<:MVector{N,T}}) where {N,T} = isbitstype(T) # see https://github.com/JuliaArrays/StaticArrays.jl/issues/27


with_eltype(::Type{<:Vector}, T) = Vector{T}
with_eltype(::Type{<:MVector{N}}, T) where {N} = MVector{N,T}
with_eltype(::Type{<:SVector{N}}, T) where {N} = SVector{N,T}
with_eltype(::Type{<:SparseVector}, T) = SparseVector{T}

# copy array and set element in one step
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




#= Miscellaneous =#

shared_sig(::OrType{<:AbstractMultivector{Sig}}...) where {Sig} = true
shared_sig(::OrType{<:AbstractMultivector}...) = false

unwrap_type(::Type{Type{T}}) where T = T
unwrap_type(T) = T

function __init__()
	if isdefined(Base.Experimental, :register_error_hint)
		Base.Experimental.register_error_hint(MethodError) do io, err, argtypes, kwargs
			# try to detect cases where the method error is due to mixing different metric signatures
			mvargtypes = unwrap_type.(filter(a -> a <: OrType{<:AbstractMultivector}, collect(argtypes)))
			if !isempty(mvargtypes) && !shared_sig(mvargtypes...)
				used_sigs = join(sprint.(show_signature, signature.(mvargtypes)), ", ", " and ")
				println(io, """
					
					Arguments have metric signatures $used_sigs.
					Perhaps the operation is not defined between these algebras?""")
			end 
		end
	end
end

const SUBSCRIPT_DIGITS   = collect("₀₁₂₃₄₅₆₇₈₉")
const SUPERSCRIPT_DIGITS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(n::Integer)   = join(SUBSCRIPT_DIGITS[begin .+ reverse(digits(n))])
superscript(n::Integer) = join(SUPERSCRIPT_DIGITS[begin .+ reverse(digits(n))])

