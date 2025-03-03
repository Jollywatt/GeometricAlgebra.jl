module NanoCAS

import Base: ==, +, -, *, /, \, ^

export WeightDict
export ProductNode, SumNode

export factor

include("weightedset.jl")
include("algebra.jl")

function variables(s::Symbol, dims::Integer...)
	indices = Iterators.product(Base.OneTo.(dims)...)
	[ProductNode(IndexNode(s, I) => 1) for I in indices]
end

struct IndexNode{N}
	symbol::Symbol
	indices::NTuple{N,Integer}
end

toexpr(a::Number) = a
toexpr(a::AbstractArray) = toexpr.(a)
toexpr(a::IndexNode) = Expr(:ref, a.symbol, a.indices...)
Base.show(io::IO, a::IndexNode) = print(io, toexpr(a))

end # module NanoCAS
