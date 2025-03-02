first_signature(::OrType{<:AbstractMultivector{Sig}}, args...) where {Sig} = Sig
first_signature(a, args...) = first_signature(args...)

"""
	use_symbolic_optim(sig) -> Bool

Whether to use symbolic optimization in algebras of metric signature `sig`.

By default, this is enabled if `dimension(sig) ≤ 8`
(in many dimensions, algebraic expressions may become too unwieldy).
"""
use_symbolic_optim(sig) = dimension(sig) <= 8

"""
	symbolic_components(label::Symbol, dims::Integer...)

Create an array of symbolic values of the specified shape.

See also [`make_symbolic`](@ref).


# Example
```jldoctest
julia> GeometricAlgebra.symbolic_components(:a, 2, 3)
2×3 Matrix{Any}:
 a[1, 1]  a[1, 2]  a[1, 3]
 a[2, 1]  a[2, 2]  a[2, 3]

julia> prod(ans)
a[1, 1]*a[1, 2]*a[1, 3]*a[2, 1]*a[2, 2]*a[2, 3]
```
"""
function symbolic_components(label::Symbol, dims::Integer...)
	var = SymbolicUtils.Sym{Array{length(dims),Real}}(label)
	indices = Iterators.product(Base.OneTo.(dims)...)
	Any[SymbolicUtils.Term{Real}(getindex, [var, I...]) for I in indices]
end

"""
	make_symbolic(a, label)

Multivector with symbolic components of the same type as the `Multivector` instance or type `a`.

See also [`symbolic_components`](@ref).

# Example

```jldoctest
julia> GeometricAlgebra.make_symbolic(Multivector{3,1}, :A)
3-component Multivector{3, 1, Vector{Any}}:
 A[1] v1
 A[2] v2
 A[3] v3
```
"""
make_symbolic(::OrType{<:Multivector{Sig,K}}, label) where {Sig,K} = Multivector{Sig,K}(symbolic_components(label, ncomponents(Multivector{Sig,K})))
make_symbolic(::OrType{F}, label) where {F<:Function} = F.instance
make_symbolic(::OrType{Val{V}}, label) where {V} = Val(V)
make_symbolic(::OrType{Type{T}}, label) where {T} = T

function toexpr(a::AbstractArray, compstype)
	SymbolicUtils.Code.toexpr(SymbolicUtils.Code.MakeArray(a, typeof(a)))
end
function toexpr(a::Multivector, compstype)
	comps = SymbolicUtils.expand.(a.comps)
	comps_expr = SymbolicUtils.Code.toexpr(SymbolicUtils.Code.MakeArray(comps, compstype))
	:( $(constructor(a))($comps_expr) )
end
toexpr(a::Tuple, compstype) = Expr(:tuple, (toexpr(ai, compstype) for ai ∈ a)...)
toexpr(a, compstype) = SymbolicUtils.Code.toexpr(a)


"""
	symbolic_multivector_eval(compstype, f, x::AbstractMultivector...)

Evaluate `f(x...)` using symbolically generated code, returning a `Multivector`
with components array of type `compstype`.

This is a generated function which first evaluates `f` on symbolic versions of
the multivector arguments `x` and then converts the symbolic result into unrolled code.

Calling `symbolic_multivector_eval(Expr, compstype, f, x...)` with `Expr` as the first argument
returns the expression to be compiled.

See also [`@symbolic_optim`](@ref).

# Example
```julia
julia> A, B = Multivector{2,1}([1, 2]), Multivector{2,1}([3, 4]);

julia> symbolic_multivector_eval(Expr, MVector, geometric_prod, A, B)
quote # prettified for readability
    let a = Multivector(args[1]).comps, b = Multivector(args[2]).comps
        comps = SymbolicUtils.Code.create_array(
            MVector, Int64, Val(1), Val((2,)),
            a[1]*b[1] + a[2]*b[2],
            a[1]*b[2] - a[2]*b[1],
        )
        Multivector{2, 0:2:2}(comps)
    end
end

julia> @btime symbolic_multivector_eval(MVector, geometric_prod, A, B);
  86.928 ns (3 allocations: 192 bytes)

julia> @btime geometric_prod(Val(:nosym), A, B); # opt-out of symbolic optim
  4.879 μs (30 allocations: 1.22 KiB)
```
"""
function symbolic_multivector_eval(::Type{Expr}, compstype::Type, f::Function, args...)
	abc = Symbol.('a' .+ (0:length(args) - 1))

	sym_args = make_symbolic.(args, abc)
	sym_result = f(sym_args...)
	I = findall(arg -> arg isa AbstractMultivector, sym_args)

	# special case for zero-component results - ensure sensible eltype
	if sym_result isa Multivector && isempty(sym_result.comps)
		T = promote_type(eltype.(args[I])...)
		return constructor(sym_result)(zeroslike(componentstype(signature(sym_result), 0, T), 0))
	end

	assignments = [:( $(abc[i]) = Multivector(args[$i]).comps ) for i in I]
	quote
		let $(assignments...)
			$(toexpr(sym_result, compstype))
		end
	end
end

@generated function symbolic_multivector_eval(compstype::Type{S}, f::Function, args...) where S
	symbolic_multivector_eval(Expr, S, f.instance, args...)
end

replace_signature(a::Multivector{Sig,K,S}, ::Val{Sig′}) where {Sig,Sig′,K,S} = Multivector{Sig′,K,S}(a.comps)
# replace_signature(a ::BasisBlade{Sig,K,T}, ::Val{Sig′}) where {Sig,Sig′,K,T} =  BasisBlade{Sig′,K,T}(a.coeff, a.bits)
replace_signature(a, ::Val) = a


function normalize_symbolic_optim_args(sig::Val, arg, args...)
	arg = arg isa BasisBlade ? Multivector(arg) : arg
	arg = replace_signature(arg, sig)
	(arg, normalize_symbolic_optim_args(sig, args...)...)
end
normalize_symbolic_optim_args(sig::Val) = ()

#=
	symbolic_optim()

Because of the rules of generated functions, we can’t call methods that may be later (re)defined
from within `symbolic_multivector_eval`. However, we still want the methods
- `dimension(sig)`
- `basis_vector_square(sig, i)`
- `componentstype(sig)`
to work for user-defined signature types, as part of the “metric signature interface”.
Since these methods may be defined in a newer world-age than `symbolic_multivector_eval`,
we must move calls to such methods outside the generated function.

To do this, the metric signature is normalized to an equivalent tuple signature, and the result of `componentstype(sig)`
is passed as an argument to — rather than being called from — `symbolic_multivector_eval`.
(We assume that `dimension(::Tuple)`, etc, are core functionality that won’t be modified by the user.)
=#

function symbolic_optim(f::Function, args...)
	# we’re replacing objects’ type parameters, so type stability is a little delicate
	Sig = first_signature(args...)
	Sig′ = canonicalize_signature(Sig)
	args′ = normalize_symbolic_optim_args(Val(Sig′), args...)
	result = symbolic_multivector_eval(componentstype(Sig, 0), f, args′...)
	replace_signature(result, Val(Sig)) # restore original signature
end

"""
	@symbolic_optim <method definition>

Convert a single method definition `f(args...)` into two methods:
1. The original method `f(Val(:nosym), args...)`, called with `Val(:nosym)` as the first argument.
   This calls the original method, without any symbolic optimization.
2. An optimized method `f(args...)` which calls [`symbolic_multivector_eval`](@ref).
   Code for this method is generated by calling `f(Val(:nosym), args...)` with symbolic versions of the `Multivector` arguments.

This is to reduce boilerplate when writing symbolically optimized versions of each method.
It only makes sense for methods with at least one `AbstractMultivector` argument for
which the exact return type is inferable.

# Example

```julia
# This macro call...
@symbolic_optim foo(a, b) = (a + b)^2
# ...is equivalent to the following two method definitions:

foo(::Val{:nosym}, a, b) = (a + b) ^ 2

function foo(a, b)
    if use_symbolic_optim(foo, a, b)
        symbolic_optim(foo, Val(:nosym), a, b)
    else
        foo(Val(:nosym), a, b)
    end
end
```
"""
macro symbolic_optim(fndef::Expr)
	fnhead, fnbody = fndef.args

	# unwrap wheres
	fnargs = fnhead::Expr
	while fnargs.head == :where
		fnargs = fnargs.args[1]::Expr
	end

	@assert fnargs.head === :call """
	`@symbolic_optim` expects a named function definition
	"""

	# give anonymous arguments (like `::Val{X}`) a name
	for arg in fnargs.args
		if arg isa Expr && arg.head === :(::) && length(arg.args) === 1
			insert!(arg.args, 1, gensym("arg"))
		end
	end

	fnhead_orig = copy(fnhead)
	fnargs_orig = copy(fnargs)
	insert!(fnargs.args, 2, :(::Val{:nosym}))

	fnname, args... = fnargs_orig.args
	quote
		$(Expr(:function, fnhead, fnbody))

		Base.@__doc__ $(Expr(:function, fnhead_orig, quote
			if use_symbolic_optim(first_signature($(args...)))
				symbolic_optim($fnname, Val(:nosym), $(args...))
			else
				$fnname(Val(:nosym), $(args...))
			end
		end))
	end |> esc
end

"""
	@symbolicga sig mv_grades expr [result_type]

Evaluate a symbolically optimised geometric algebra expression.

Upon macro expansion, `expr` is evaluated with symbolic multivectors (specified by `mv_grades`)
in the algebra defined by metric signature `sig`. The resulting symbolic expression
is then compiled and executed at runtime.

The `mv_grades` argument is a `NamedTuple` where `keys(mv_grades)`
defines the identifiers in `expr` to be interpreted as `Multivector`s,
while `values(mv_grades)` defines their respective grades.
The identifiers must exist at runtime, and can be a `Multivector` with matching
signature/grade or any iterable with the correct number of components.

If `result_type` is given, then the components of the resulting multivector
are converted to that type. The result type `T` should implement `T(::Tuple)`,
e.g., `Tuple` or `MVector`.

!!! warning
	Operations that are not amenable to symbolic evaluation
	(such as `exp`, `log`, `sqrt`, etc) are not supported.

	(You may test if operations work on symbolic multivectors
	created with [`GeometricAlgebra.make_symbolic`](@ref).)

# Examples
```jldoctest
julia> v = (1, 2, 0); R = exp(Multivector{3,2}([0, π/8, 0]));

julia> # Rotate a tuple (interpreted as a grade 1 vector)
       # by a rotor, returning a tuple.
       @symbolicga 3 (v=1, R=0:2:4) grade(R*v*~R, 1) Tuple
(0.7071067811865475, 2.0, -0.7071067811865476)
```
```julia
# This macro call...
@symbolicga 3 (a=1, b=1) a∧b
# ...is equivalent to the following:
let a = Multivector{3, 1}(a).comps, b = Multivector{3, 1}(b).comps
    Multivector{3, 2}(MVector(
        a[1]*b[2] - a[2]*b[1],
        a[1]*b[3] - a[3]*b[1],
        a[2]*b[3] - a[3]*b[2],
    ))
end
```
"""
macro symbolicga(sig, mv_grades, expr, result_type=nothing)
	Sig = eval(sig)
	mv_grades = eval(mv_grades)
	@assert mv_grades isa NamedTuple
	result_type = eval(result_type)

	labels = keys(mv_grades)

	symbolic_assignments = map(labels, mv_grades) do label, K
		mv = make_symbolic(Multivector{Sig,K}, label)
		:($label = $mv)
	end
	symbolic_result = @eval let $(symbolic_assignments...)
		$expr
	end

	if isnothing(result_type)
		result_expr = toexpr(symbolic_result, componentstype(Sig, 0))
	else
		symbolic_result isa Multivector || error("@symbolicga expression must evaluate to a Multivector when result_type ($result_type) is given; got $(typeof(symbolic_result)).")
		symbolic_comps = Tuple(toexpr(comp, nothing) for comp in symbolic_result.comps)
		result_expr = :($result_type(($(symbolic_comps...),)))
	end

	expr_assignments = map(labels, mv_grades) do label, K
		:($label = Multivector{$Sig,$K}($label).comps)
	end

	quote
		let $(expr_assignments...)
			$result_expr
		end
	end |> esc

end
