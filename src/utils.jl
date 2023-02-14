#= Tools to handle symbolic eltypes =#

#=
Multivectors with purely symbolic components have SymbolicUtils.Symbolic eltype.
For sym::SymbolicUtils.Sym, iszero(sym) is a symbolic expression “sym == 0”,
but we generally want `iszero(sym) == false` because it *could* be non-zero.
For this purpose, we have `isnumberzero` and friends.
Furthermore, we want to support eltypes of `Any` (e.g., when combining a symbolic
multivector with eltype `Int`), so need `numberzero(::Any)` to work.
=#


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

# initialize zeros while abstracting away array type
zeroslike(::Type{<:AbstractArray{T}}, dims...) where {T<:Number} = zeros(T, dims...)
zeroslike(::Type{<:AbstractArray{T}}, dims...) where {T} = convert(Array{Any}, zeros(Int, dims...))
zeroslike(::Type{<:MArray{N,T}}, dims...) where {N,T<:Number} = zeros(MArray{Tuple{dims...},T})
zeroslike(::Type{<:MArray{N,T}}, dims...) where {N,T} = MArray{Tuple{dims...},Any}(zeros(MArray{Tuple{dims...},Int}))
zeroslike(::Type{<:SArray{N,T}}, dims...) where {N,T<:Number} = zeros(SArray{Tuple{dims...},T})
zeroslike(::Type{<:SArray{N,T}}, dims...) where {N,T} = SArray{Tuple{dims...},Any}(zeros(SArray{Tuple{dims...},Int}))
zeroslike(::Type{<:SparseVector{Tv}}, dims...) where {Tv} = spzeros(Tv, dims...)
zeroslike(::Type{<:SubArray{T,N,P}}, dims...) where {T,N,P} = zeroslike(P, dims...)


@static if VERSION < v"1.7"
	ismutabletype(a::DataType) = a.mutable
end

# whether array supports `setindex!`
issetindexable(::Type{<:AbstractSparseArray}) = true
issetindexable(::Type{<:MVector{N,T}}) where {N,T} = isbitstype(T) # see https://github.com/JuliaArrays/StaticArrays.jl/issues/27
issetindexable(T::Type) = ismutabletype(T)
issetindexable(a) = issetindexable(typeof(a))




#= Miscellaneous =#

shared_sig(::OrType{<:AbstractMultivector{Sig}}...) where {Sig} = true
shared_sig(::OrType{<:AbstractMultivector}...) = false

function __init__()
	if isdefined(Base.Experimental, :register_error_hint)
		Base.Experimental.register_error_hint(MethodError) do io, err, argtypes, kwargs
			# try to detect cases where the method error is due to mixing different metric signatures
			mvargtypes = filter(a -> a <: AbstractMultivector, collect(argtypes))
			if !isempty(mvargtypes) && !shared_sig(mvargtypes...)
				used_sigs = join(sprint.(show_signature, signature.(mvargtypes)), ", ", " and ")
				println(io, """
					
					Arguments have metric signatures $used_sigs.
					Perhaps the operation is not defined between these algebras?""")
			end 
		end
	end

	# init_symbolic_optim()
end

const SUBSCRIPT_DIGITS   = collect("₀₁₂₃₄₅₆₇₈₉")
const SUPERSCRIPT_DIGITS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(n::Integer)   = join(SUBSCRIPT_DIGITS[begin .+ reverse(digits(n))])
superscript(n::Integer) = join(SUPERSCRIPT_DIGITS[begin .+ reverse(digits(n))])

