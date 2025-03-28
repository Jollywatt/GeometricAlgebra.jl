struct SubexprPointer
	x::UInt
end

"""
	SubexprList <: AbstractDict{SubexprPointer,Expr}

Represents a directed acyclic graph of expressions as an ordered list of definitions.
The last entry defines the full expression in terms of preceding subexpressions.
Subexpressions are referenced with a `SubexprPointer`, which are the keys of the underlying `OrderedDict`
and are pretty printed as Greek letters.

See also [`subexprs`](@ref).

# Examples
```jldoctest
julia> subexprs(:(A + f(A) + g(f(A))))
SubexprList with 3 entries:
  α => :(f(A))
  β => :(g(α))
  γ => :(A + α + β)
```
"""
struct SubexprList <: AbstractDict{SubexprPointer,Expr}
	defs::OrderedDict{SubexprPointer,Expr}
end
SubexprList() = SubexprList(OrderedDict{SubexprPointer,Expr}())

Base.iterate(a::SubexprList, args...) = iterate(a.defs, args...)
Base.length(a::SubexprList) = length(a.defs)
Base.last(a::SubexprList) = last(a.defs)
Base.setindex!(a::SubexprList, args...) = setindex!(a.defs, args...)
Base.get(a::SubexprList, k, default) = get(a.defs, k, default)


function letter(i::Integer)
	alphabet = collect("αβγδεζηθκλμνξπρςστυφχψωΓΔΘΛΞΠΣΦΨΩ")
	n = length(alphabet) 
	base = alphabet[mod1(i, n)]
	Symbol(base^cld(i, n))
end

function Base.show(io::IO, ::MIME"text/plain", s::SubexprList)
	names = Dict(k => letter(i) for (i, k) in enumerate(keys(s.defs)))
	io = IOContext(io, :pointernames => names)
	@invoke show(io, MIME("text/plain"), s::AbstractDict)
end

function Base.show(io::IO, s::SubexprPointer)
	if :pointernames in keys(io)
		printstyled(io, io[:pointernames][s], bold=true)
	else
		@invoke show(io, s::Any)
	end
end

subexprs!(l::SubexprList, a::Any) = a
function subexprs!(l::SubexprList, expr::Expr)
	ref = SubexprPointer(hash(expr))
	ref in keys(l) && return ref

	args = map(expr.args) do arg
		subexprs!(l, arg)
	end

	expr.head == :ref && return expr

	new = Expr(expr.head, args...)
	push!(l, ref => new)
	ref
end

"""
	subexprs(::Expr)::SubexprList

Flatten an expression tree into a list of atomic expressions so that common subexpressions are identified.

See also [`squash`](@ref).

# Example
```jldoctest
julia> subexprs(:(A + f(A) + g(f(A))))
SubexprList with 3 entries:
  α => :(f(A))
  β => :(g(α))
  γ => :(A + α + β)

julia> toexpr(ans, pretty=true)
:(let α = f(A), β = g(α)
      A + α + β
  end)
```
"""
function subexprs(expr::Expr)
	l = SubexprList()
	subexprs!(l, expr)
	l
end

countrefs!(counts, ::Any) = nothing
function countrefs!(counts, ref::SubexprPointer)
	counts[ref] = get(counts, ref, 0) + 1
end
function countrefs!(counts, expr::Expr)
	for arg in expr.args
		countrefs!(counts, arg)
	end
end
function countrefs!(counts, l::SubexprList)
	for (ref, v) in l
		countrefs!(counts, v)
	end
end

function countrefs(a)
	counts = Dict{SubexprPointer,Int}()
	countrefs!(counts, a)
	counts
end


substitute(a, ::Any) = a
substitute(ref::SubexprPointer, subs::Dict{SubexprPointer}) = get(subs, ref, ref)
function substitute(expr::Expr, subs)
	args = map(expr.args) do arg
		substitute(arg, subs)
	end
	Expr(expr.head, args...)
end


"""
	squash(::SubexprList, maxcount=1)::SubexprList

Eliminate subexpressions which are only referenced up to `maxcount` times
substituting their definitions into subsequent expressions.

# Example

```jldoctest
julia> subexprs(:(A + f(A) + g(f(A))^2))
SubexprList with 4 entries:
  α => :(f(A))
  β => :(g(α))
  γ => :(β ^ 2)
  δ => :(A + α + γ)

julia> squash(ans, 1)
ERROR: UndefVarError: `squash` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1

julia> squash(ans, 2)
ERROR: UndefVarError: `squash` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Stacktrace:
 [1] top-level scope
   @ none:1
```
"""
function squash(l::SubexprList, maxcount=1)
	counts = countrefs(l)
	tosquash = keys(filter(<=(maxcount)∘last, counts))
	subs = filter(in(tosquash)∘first, l)
	out = SubexprList()
	for (ref, v) in l
		new = substitute(v, subs)
		if ref ∈ keys(subs)
			subs[ref] = new
		end
		if ref ∉ tosquash
			out[ref] = new
		end
	end
	out
end


"""
	cse(a)

Perform [common subexpression elimination](https://en.wikipedia.org/wiki/Common_subexpression_elimination)
on an `Expr` or value (to be converted with [`toexpr`](@ref)).

Returns a `let ... end` block.

```jldoctest
julia> cse(:(A + f(A) + g(f(A))^2))
:(let α = f(A)
      A + α + g(α) ^ 2
  end)
```
"""
cse = toexpr∘squash∘subexprs∘toexpr
