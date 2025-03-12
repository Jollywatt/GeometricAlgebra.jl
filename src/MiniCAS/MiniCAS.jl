module MiniCAS

import Base: ==, +, -, *, /, \, ^
using OrderedCollections: OrderedDict

export ProductNode, SumNode
export SubexprList

export variable, variables, factor, toexpr, subexprs, cse

toexpr(a::Union{Symbol,Expr,Number}) = a

include("weightedset.jl")
include("algebra.jl")
include("cse.jl")
include("toexpr.jl")


variable(label::Symbol) = ProductNode(label => 1)
function variables(s::Symbol, dims::Integer...)
	indices = Iterators.product(Base.OneTo.(dims)...)
	[ProductNode(Expr(:ref, s, I...) => 1) for I in indices]
end

"""
	variable(label::Symbol)
	variables(label::Symbol, dims::Integer...)

A symbolic real value or array of symbolic components.


# Example
```jldoctest
julia> variables(:a, 2, 2)
2×2 Matrix{ProductNode{Expr}}:
 a[1, 1]  a[1, 2]
 a[2, 1]  a[2, 2]

julia> prod(ans) + variable(:b)
SumNode{Any, Int64}:
 a[1, 1] * a[1, 2] * a[2, 1] * a[2, 2] + b
```
"""
variable, variables



end # module MiniCAS
