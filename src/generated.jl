
"""
	use_symbolic_optim(sig) -> Bool

Whether to use symbolic optimization in algebras of metric signature `sig`.

By default, this is enabled if `dimension(sig) ≤ 8`
(in many dimensions, algebraic expressions may become too unwieldy).
"""
use_symbolic_optim(sig) = dimension(sig) <= 8
use_symbolic_optim(::Val{Sig}) where Sig = use_symbolic_optim(Sig)

"""
	Multivector{Sig,K}(sym::Symbol)

Multivector with independent symbolic components.

See also [`make_symbolic`](@ref).

# Example
```jldoctest
julia> a = Multivector{3,1}(:a)
3-component Multivector{3, 1, Vector{ProductNode{Expr}}}:
 a[1] v1
 a[2] v2
 a[3] v3

julia> a ⊙ a
SumNode{Expr, Int64}:
 a[1] ^ 2 + a[2] ^ 2 + a[3] ^ 2
```
"""
Multivector{Sig,K}(sym::Symbol) where {Sig,K} = make_symbolic(Multivector{Sig,K}, sym)

"""
	make_symbolic(a, label)

Multivector with symbolic components of the same type as the `Multivector` instance or type `a`.

See also [`MiniCAS.variables`](@ref).

# Example

```jldoctest
julia> GeometricAlgebra.make_symbolic(Multivector{3,1}, :A)
3-component Multivector{3, 1, Vector{ProductNode{Expr}}}:
 A[1] v1
 A[2] v2
 A[3] v3

```
"""
make_symbolic(::OrType{<:Multivector{Sig,K}}, label) where {Sig,K} = Multivector{Sig,K}(MiniCAS.variables(label, ncomponents(Multivector{Sig,K})))
make_symbolic(::OrType{F}, label) where {F<:Function} = F.instance
make_symbolic(::OrType{Val{V}}, label) where {V} = Val(V)
# make_symbolic(::OrType{Type{T}}, label) where {T} = T
make_symbolic(::OrType{<:Number}, label) = MiniCAS.variable(label)


import .MiniCAS: toexpr, factor

toexpr(a::Multivector{Sig}) where Sig = toexpr(a, Val(Sig))

toexpr(a::Multivector, ::Val{Sig}) where Sig = :(Multivector{$Sig,$(grade(a))}($(toexpr.(a.comps)...)))
toexpr(a, ::Val) = toexpr(a)

toexpr(a::Vector) = :([$(a...)])
toexpr(a::Matrix) = :(Base.hvcat($(size(a, 2)), $(toexpr.(permutedims(a))...)))
toexpr(a::SArray{Size}) where Size = :(SArray{$Size}($(toexpr.(a)...)))
toexpr(a::MArray{Size}) where Size = :(MArray{$Size}($(toexpr.(a)...)))

factor(a::Multivector) = constructor(a)(factor.(a.comps))
factor(a::Number) = a


"""
	symbolic_multivector_eval(sig::Val, f::Function, args...)

Evaluate `f(args...)` using symbolically optimised code for operations on `Multivector`s.

This is a generated function which first evaluates `f` on symbolic versions of
the multivector arguments `make_symbolic.(args)` and then converts the symbolic result into unrolled code.

If the result is a `Multivector`, it is given the metric signature `sig`.

Calling `symbolic_multivector_eval(Expr, sig, f, args...)` with `Expr` as the first argument
returns the expression to be compiled.

See also [`@symbolic_optim`](@ref).

# Example
```julia
julia> A, B = randn(Multivector{3,0:3}, 2)

julia> symbolic_multivector_eval(Expr, Val(2), geometric_prod, A, B)
:(let a = (Multivector(args[1])).comps, b = (Multivector(args[2])).comps
      Multivector{2, 0:2}(
          a[1] * b[1] + -1 * (b[4] * a[4]) + b[2] * a[2] + b[3] * a[3],
          b[3] * a[4] + b[1] * a[2] + -1 * (b[4] * a[3]) + a[1] * b[2],
          -1 * (a[4] * b[2]) + b[4] * a[2] + b[3] * a[1] + b[1] * a[3],
          b[3] * a[2] + a[4] * b[1] + -1 * (b[2] * a[3]) + b[4] * a[1],
      )
  end)

julia> @btime symbolic_multivector_eval(Val(2), geometric_prod, A, B);
  19.684 ns (2 allocations: 64 bytes)

julia> @btime geometric_prod(Val(:nosym), A, B); # opt-out of symbolic optim
  7.312 μs (125 allocations: 4.67 KiB)
```
"""
function symbolic_multivector_eval(::Type{Expr}, sig::Val, f::Function, args...; simplify=true)
	abc = Symbol.('a' .+ (0:length(args) - 1))

	sym_args = make_symbolic.(args, abc)
	sym_result = f(sym_args...)
	I = findall(arg -> arg isa AbstractMultivector, sym_args)

	if sym_result isa Multivector
		# special case for zero-component results - ensure sensible eltype
		if isempty(sym_result.comps)
			T = promote_type(eltype.(args[I])...)
			comps = zeroslike(componentstype(signature(sym_result), 0, T), 0)
			return constructor(sym_result)(comps)
		end

		sym_result = MiniCAS.factor(sym_result)
	end

	expr = toexpr(sym_result, sig)
	if simplify && expr isa Expr
		expr = MiniCAS.cse(expr)
	end

	assignments = [i in I ? :( $(abc[i]) = Multivector(args[$i]).comps ) : :($(abc[i]) = args[$i]) for i in eachindex(args)]
	:(let $(assignments...)
		$expr
	end)
end


@generated function symbolic_multivector_eval(::Val{Sig}, f::Function, args...) where Sig
	@assert isdefined(f, :instance)
	1
	symbolic_multivector_eval(Expr, Val(Sig), f.instance, args...)
end


# return the first metric signature parameter encountered in an argument list
first_signature(::OrType{<:AbstractMultivector{Sig}}, args...) where Sig = Val(Sig)
first_signature(a, args...) = first_signature(args...)

replace_signature(a::Multivector{Sig,K,S}, ::Val{Sig′}) where {Sig,Sig′,K,S} = Multivector{Sig′,K,S}(a.comps)
replace_signature(a, ::Val) = a

canonicalize(a::Multivector) = replace_signature(a, Val(canonical_signature(signature(a))))
canonicalize(a::BasisBlade) = canonicalize(Multivector(a))
canonicalize(a) = a

"""
	symbolic_optim(f, args...)

Evaluate `f(args...)` by invoking the generated function [`symbolic_multivector_eval`](@ref).

Because of the rules of generated functions, [`symbolic_multivector_eval`](@ref) must not call
methods that may be later (re)defined. However, we still want the methods
- `dimension(sig)`
- `basis_vector_square(sig, i)`
- `componentstype(sig)`
to work for user-defined signature types, as part of the “metric signature interface”.
Since these methods may be defined in a newer world-age than `symbolic_multivector_eval`,
we must move calls to such methods outside the generated function.

To do this, the metric signatures in `args` are replaced with the equivalent canonical tuple signature.
(We assume that `dimension(::Tuple)`, etc, are core functionality that won’t be modified by the user.)

!!! warning
	If `f(args...)` is a `Multivector`, its signature is assumed to be identical to the signature
	of the first `AbstractMultivector` argument in `args`.
	(The actual signature is lost because signatures are converted to canonical tuples.)
"""
function symbolic_optim(f::Function, args...)
	# we’re replacing objects’ type parameters, so type stability is a little delicate
	args′ = map(canonicalize, args)
	Sig::Val = first_signature(args...) # guess the original (non-canonical) signature of f(args...)
	result = symbolic_multivector_eval(Sig, f, args′...)
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
@symbolicga 3 (a=1, b=1) wedge(a, b)
# ...is equivalent to the following:
let a = Multivector{3, 1}(a).comps, b = Multivector{3, 1}(b).comps
    Multivector{3, 2}(
        a[1]*b[2] - a[2]*b[1],
        a[1]*b[3] - a[3]*b[1],
        a[2]*b[3] - a[3]*b[2],
    )
end
```
"""
macro symbolicga(sig, mv_grades, expr, result_type=nothing)
	Sig = eval(sig)
	mv_grades = eval(mv_grades)
	@assert mv_grades isa NamedTuple
	labels = keys(mv_grades)
	result_type = eval(result_type)


	symbolic_assignments = map(labels, mv_grades) do label, K
		mv = make_symbolic(Multivector{Sig,K}, label)
		:($label = $mv)
	end
	symbolic_result = @eval let $(symbolic_assignments...)
		$expr
	end

	if isnothing(result_type)
		result_expr = toexpr(symbolic_result)
	else
		symbolic_result isa Multivector || throw(ArgumentError("""
			@symbolicga expression must evaluate to a Multivector when result_type=$result_type is specified; \
			got a value of type $(typeof(symbolic_result))."""))

		result_expr = :($makevec($result_type, $(toexpr.(symbolic_result.comps)...)))

	end

	result_expr = MiniCAS.cse(result_expr)

	comps_assignments = map(labels, mv_grades) do label, K
		:($label = Multivector{$Sig,$K}($label).comps)
	end

	quote
		let $(comps_assignments...)
			$result_expr
		end
	end |> esc

end
