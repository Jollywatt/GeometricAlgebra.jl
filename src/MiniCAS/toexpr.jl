function toexpr(a::ProductNode)
	terms = [isone(v) ? toexpr(k) : Expr(:call, :^, toexpr(k), v) for (k, v) in a.x]
	isempty(terms) && return 1
	isone(length(terms)) && return only(terms)
	Expr(:call, :*, terms...)
end

function toexpr(a::SumNode)
	function pretty(v, expr)
		v == 1 && return expr
		v == -1 && return :(-$expr)
		:($v*$expr)
	end
	terms = [pretty(v, toexpr(k)) for (k, v) in a.x]
	isempty(terms) && return 0
	isone(length(terms)) && return only(terms)
	Expr(:call, :+, terms...)
end


#= show methods =#

Base.show(io::IO, a::Union{Π,Σ}) = print(io, toexpr(a))
Base.show(io::IO, ::MIME"text/plain", a::Union{Π,Σ}) = print(io, typeof(a), ":\n ", a)

debug(a) = a
debug(a::Π) = Expr(:call, :Π, (:($(debug(k))^$v) for (k, v) in a.x)...)
debug(a::Σ) = Expr(:call, :Σ, (:($v*$(debug(k))) for (k, v) in a.x)...)



struct IndexNode{N}
	symbol::Symbol
	indices::NTuple{N,Integer}
end
toexpr(a::IndexNode) = Expr(:ref, a.symbol, a.indices...)
Base.show(io::IO, a::IndexNode) = print(io, toexpr(a))


"""
	toexpr(::SubexprList; pretty=false)

Render a subexpression list as a `let ... end` expression.

# Example
```jldoctest
julia> toexpr(subexprs(:(x^2 + f(x^2))), pretty=true)
:(let α = x ^ 2, β = f(α)
      α + β
  end)
```
"""
function toexpr(l::SubexprList; pretty=true)
	names = Dict(k => pretty ? letter(i) : gensym() for (i, k) in enumerate(keys(l.defs)))
	c = collect(l)
	defs = [names[k] => substitute(v, names) for (k, v) in l]
	lets = [:($k = $v) for (k, v) in defs[1:end-1]]
	Expr(:let, Expr(:block, lets...), last(defs[end]))
end
