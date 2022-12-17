"""
	use_symbolic_optim(sig) -> Bool

Whether to use symbolic code generation to optimize operations
in algebras of metric signature `sig`.

By default, this is enabled if `dimension(sig) ≤ 8` as a heuristic
(in many dimensions, algebraic expressions may become too unwieldy).
"""
use_symbolic_optim(sig) = dimension(sig) <= 8
use_symbolic_optim(::Function, ::Union{Scalar,Function,AbstractMultivector{Sig}}...) where {Sig} = use_symbolic_optim(Sig)

function symbolic_components(label, dims...)
	var = SymbolicUtils.Sym{Array{length(dims),Real}}(label)
	indices = Iterators.product(Base.OneTo.(dims)...)
	Any[SymbolicUtils.Term{Real}(getindex, [var, I...]) for I in indices]
end

symbolic_argument(::OrType{<:BasisBlade{Sig,K}}, label) where {Sig,K} = symbolic_argument(Multivector{Sig,K}, label)
symbolic_argument(A::OrType{<:Multivector{Sig,K}}, label) where {Sig,K} = Multivector{Sig,K}(symbolic_components(label, ncomponents(A)))
symbolic_argument(::OrType{F}, label) where {F<:Function} = F.instance
symbolic_argument(::OrType{Val{V}}, label) where {V} = Val(V)

toexpr(a, T) = SymbolicUtils.Code.toexpr(a)
function toexpr(a::AbstractMultivector, T)
	compstype = componentstype(signature(a), ncomponents(a), T)
	comps_expr = SymbolicUtils.Code.toexpr(SymbolicUtils.Code.MakeArray(a.comps, compstype, T))

	:( $(constructor(a))($comps_expr) )
end


"""
	symbolic_multivector_eval(f, x::AbstractMultivector...)

Evaluate `f(x...)` using symbolically generated code.

This is a generated function which first evaluates `f` on symbolic versions of
the multivector arguments `x` and then converts the symbolic result into unrolled code.

Calling `symbolic_multivector_eval(Expr, f, x...)` with `Expr` as the first argument
returns the unevaluated code as an expression (for introspection).
"""
@generated function symbolic_multivector_eval(f::Function, args...)
	symbolic_multivector_eval(Expr, f.instance, args...)
end
function symbolic_multivector_eval(::Type{Expr}, f::Function, args...)
	abc = Symbol.('a' .+ (0:length(args) - 1))

	sym_args = symbolic_argument.(args, abc)
	sym_result = f(sym_args...)

	I = findall(T -> T isa Multivector, sym_args)
	assignments = [:( $(abc[i]) = components(args[$i]) ) for i in I]
	T = numberorany(promote_type(eltype.(args[I])...))

	quote
		let $(assignments...)
			$(toexpr(sym_result, T))
		end
	end
end



"""
	@symbolic_optim

Applied to a method definition accepting `AbstractMultivector` arguments,
define an optimized method which calls `symbolic_multivector_eval`
in addition to the original method.

The original, unoptimized method is called with `Val(:nosym)` as the first argument.

# Example

```julia
@symbolic_optim foo(a, b) = (a + b)^2
# ...is equivalent to defining the two methods below:

foo(::Val{:nosym}, a, b) = (a + b) ^ 2

function foo(a, b)
    if use_symbolic_optim(foo, a, b)
        symbolic_multivector_eval(foo, Val(:nosym), a, b)
    else
        foo(Val(:nosym), a, b)
    end
end
```
"""
macro symbolic_optim(fndef::Expr)
	fnhead, fnbody = fndef.args

	fnargs = fnhead::Expr
	while fnargs.head == :where
		fnargs = fnargs.args[1]::Expr
	end
	fnhead_orig = copy(fnhead)
	fnargs_orig = copy(fnargs)
	insert!(fnargs.args, 2, :(::Val{:nosym}))

	fnname, mvargs... = fnargs_orig.args
	quote
		$(Expr(:function, fnhead, fnbody))

		Base.@__doc__ $(Expr(:function, fnhead_orig, quote
			if use_symbolic_optim($fnname, $(mvargs...))
				symbolic_multivector_eval($fnname, Val(:nosym), $(mvargs...))
			else
				$fnname(Val(:nosym), $(mvargs...))
			end
		end))
	end |> esc
end


# way to convert a BasisBlade to a Multivector without allocating a full components array
# TODO: take this more seriously
function components(a::BasisBlade{Sig,K}) where {Sig,K}
	i = findfirst(==(bitsof(a)), componentbits(Val(dimension(Sig)), Val(K)))
	SingletonVector(a.coeff, i, ncomponents(Sig, K))
end
components(a::Multivector) = a.comps

struct SingletonVector{T} <: AbstractVector{T}
	el::T
	index::Int
	length::Int
end
Base.length(a::SingletonVector) = a.length
Base.size(a::SingletonVector) = (length(a),)
Base.getindex(a::SingletonVector{T}, i) where {T} = a.index == i ? a.el : numberzero(T)
function Base.iterate(a::SingletonVector, i = 1)
	i > a.length && return
	(a[i], i + 1)
end
