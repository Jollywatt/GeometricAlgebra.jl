#= Tools to handle symbolic eltypes =#

#=
Multivectors with purely symbolic components have SymbolicUtils.Symbolic eltype.
For sym::SymbolicUtils.Sym, iszero(sym) is a symbolic expression “sym == 0”,
but we generally want `iszero(sym) == false` because it *could* be non-zero.
For this purpose, we have `isnumberzero` and friends.
Furthermore, we want to support eltypes of `Any` (e.g., when combining a symbolic
multivector with eltype `Int`), so need `numberzero(::Any)` to work.
=#


isnumberzero(a) = iszero(a)
isnumberone(a) = isone(a)

# don't treat a symbolic type as one/zero, even if it could represent one
# isnumberzero(a) = false
# isnumberone(a) = false

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
zeroslike(a::AbstractArray, dims...) = zeroslike(typeof(a), dims...)


@static if VERSION < v"1.7"
	ismutabletype(a::DataType) = a.mutable
end

# whether array supports `setindex!`
issetindexable(::Type{<:AbstractSparseArray}) = true
issetindexable(::Type{<:MVector{N,T}}) where {N,T} = isbitstype(T) # see https://github.com/JuliaArrays/StaticArrays.jl/issues/27
issetindexable(T::Type) = ismutabletype(T)
issetindexable(a) = issetindexable(typeof(a))


makevec(::Type{<:Tuple}, comps...) = comps
makevec(::Type{<:Vector}, comps...) = collect(comps)
makevec(::Type{<:MVector}, comps...) = MVector(comps...)
makevec(::Type{<:SVector}, comps...) = SVector(comps...)

"""
	SingletonVector(el, index, length)

Efficient representation of a vector of all zeros
except for the single element `el` at the given index.
"""
struct SingletonVector{T} <: AbstractVector{T}
	el::T
	index::Int
	length::Int
end

Base.length(a::SingletonVector) = a.length
Base.size(a::SingletonVector) = (length(a),)
Base.eltype(::SingletonVector{T}) where {T} = numberorany(T)
Base.iszero(a::SingletonVector) = iszero(a.el)

Base.getindex(a::SingletonVector{T}, i::Integer) where {T} = a.index == i ? a.el : numberzero(T)
Base.getindex(a::SingletonVector, I::UnitRange) = SingletonVector(a.el, a.index - first(I) + 1, length(I))

function Base.iterate(a::SingletonVector, i = 1)
	i > a.length && return
	(a[i], i + 1)
end
Base.:(==)(a::SingletonVector, b::AbstractVector) = length(a) == length(b) && mapreduce(==, &, a, b)

for op in [:*, :/]
	@eval Base.$op(a::SingletonVector, b::Number) = SingletonVector($op(a.el, b), a.index, a.length)
	@eval Base.$op(a::Number, b::SingletonVector) = SingletonVector($op(a, b.el), b.index, b.length)
end
Base.:(//)(a::SingletonVector, b::Number) = SingletonVector(a.el//b, a.index, a.length)


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
end

const SUBSCRIPT_DIGITS   = collect("₀₁₂₃₄₅₆₇₈₉")
const SUPERSCRIPT_DIGITS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(n::Integer)   = (n < 0 ? '₋' : "")*join(SUBSCRIPT_DIGITS[begin .+ reverse(digits(abs(n)))])
superscript(n::Integer) = (n < 0 ? '⁻' : "")*join(SUPERSCRIPT_DIGITS[begin .+ reverse(digits(abs(n)))])

