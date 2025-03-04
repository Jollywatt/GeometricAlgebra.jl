module MiniCAS

import Base: ==, +, -, *, /, \, ^
using OrderedCollections: OrderedDict

export WeightDict
export ProductNode, SumNode

export variable, variables, factor, toexpr

include("weightedset.jl")
include("algebra.jl")
include("cse.jl")

variable(s::Symbol) = ProductNode(s => 1)

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

end # module MiniCAS
