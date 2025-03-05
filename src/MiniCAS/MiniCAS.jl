module MiniCAS

import Base: ==, +, -, *, /, \, ^
using OrderedCollections: OrderedDict

export WeightDict
export ProductNode, SumNode

export variable, variables, factor, toexpr, subexprs, cse

toexpr(a::Union{Symbol,Expr,Number}) = a

include("weightedset.jl")
include("algebra.jl")
include("cse.jl")
include("toexpr.jl")

variable(s::Symbol) = ProductNode(s => 1)

function variables(s::Symbol, dims::Integer...)
	indices = Iterators.product(Base.OneTo.(dims)...)
	[ProductNode(IndexNode(s, I) => 1) for I in indices]
end




end # module MiniCAS
