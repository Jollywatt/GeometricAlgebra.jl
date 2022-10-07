# fallback herustic: use symbolic optimisation unless dim too big
use_symbolic_optim(sig) = dimension(sig) <= 8


function symbolic_components(label, dims...)
	var = SymbolicUtils.Sym{Array{length(dims),Real}}(label)
	indices = Iterators.product(Base.OneTo.(dims)...)
	[SymbolicUtils.Term{Real}(getindex, [var, I...]) for I in indices]
end


symbolic_multivector(a::Type{<:Blade{Sig,K}}, label) where {Sig,K} = symbolic_multivector(Multivector{Sig,K}, label)
function symbolic_multivector(A::Type{<:CompositeMultivector{Sig,C}}, label) where {Sig,C}
	constructor(A)(symbolic_components(label, ncomponents(A)))
end
symbolic_multivector(A::Type{Blade{Sig,K,T}}, label) where {Sig,K,T} = symbolic_multivector(Multivector{Sig,K,componentstype(Sig, ncomponents(Sig, K), T)}, label)
symbolic_multivector(a::AbstractMultivector, label) = symbolic_multivector(typeof(a), label)


"""
	symbolic_optim(f, a, b, ...)

Trace evaluation of `f(a, b, ...)::CompositeMultivector` on symbolic versions of each
`AbstractMultivector` instance or type `a`, `b`, ..., returning an expression suitable
as the body of `@generated` function.

!!! important
	The names of the function arguments must be `"a"`, `"b"`, `"c"`, etc, as these are
	the names used in the expression retuned by `symbolic_optim`.

If `use_symbolic_optim(sig)` returns `false`, the function body simply calls `f(a, b, ...)`.

# Examples
```jldoctest
using MacroTools: prettify
u, v = Multivector.(basis(2))
ex = GeometricAlgebra.symbolic_optim(*, u, v) |> prettify

# output
:(let a = components(a), b = components(b)
      comps = create_array(Vector{Any}, Int64, Val{1}(), Val{(4,)}(), a[1] * b[1] + a[2] * b[2], 0, 0, a[1] * b[2] + (-1 * a[2]) * b[1])
      (MixedMultivector{2})(comps)
  end)
```
"""
function symbolic_optim(f, x::OrType{<:AbstractMultivector{Sig}}...) where {Sig}
	abc = Symbol.('a' .+ (0:length(x) - 1))

	if !use_symbolic_optim(Sig)
		return :( $f($(abc...)) )
	end
	
	x_symb = symbolic_multivector.(x, abc)
	y_symb = f(x_symb...)
	comps = y_symb.comps

	T = numberorany(promote_type(eltype.(x)...))
	comps_expr = SymbolicUtils.Code.MakeArray(comps, typeof(comps), T)

	assignments = [:( $a = components($a) ) for a in abc]
	quote
		let $(assignments...)
			comps = $(SymbolicUtils.Code.toexpr(comps_expr))
			$(constructor(y_symb))(comps)
		end
	end
end

# way to convert a Blade to a Multivector without allocating a full components array
# TODO: take this more seriously
function components(a::Blade{Sig,K}) where {Sig,K}
	i = bits_to_mv_index(bitsof(a))
	SingletonVector(a.coeff, i, ncomponents(Sig, K))
end
components(a::CompositeMultivector) = a.comps

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
