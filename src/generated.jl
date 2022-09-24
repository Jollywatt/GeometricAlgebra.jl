function symbolic_components(label, dims...)
	var = SymbolicUtils.Sym{Array{length(dims),Real}}(label)
	indices = Iterators.product(Base.OneTo.(dims)...)
	[SymbolicUtils.Term{Real}(getindex, [var, I...]) for I in indices]
end



with_eltype(::Type{<:Multivector{Sig,K,C}}, T) where {Sig,K,C} = Multivector{Sig,K,with_eltype(C, T)}
with_eltype(::Type{<:MixedMultivector{Sig,C}}, T) where {Sig,C} = MixedMultivector{Sig,with_eltype(C, T)}

function symbolic_multivector(A::Type{<:CompositeMultivector{Sig,C}}, label) where {Sig,C}
	constructor(A)(symbolic_components(label, ncomponents(A)))
end
symbolic_multivector(A::Type{Blade{Sig,K,T}}, label) where {Sig,K,T} = symbolic_multivector(Multivector{Sig,K,componentstype(Sig, ncomponents(Sig, K), T)}, label)
symbolic_multivector(a::AbstractMultivector, label) = symbolic_multivector(typeof(a), label)


"""
	generated_multivector_function(f, a, b, ...)

Trace evaluation of `f(a, b, ...)::CompositeMultivector` on symbolic versions of each
`AbstractMultivector` instance or type `a`, `b`, ..., returning an expression body
suitable for a `@generated` function.

The names of the arguments to be passed to the function
body are the literal symbols `:a`, `:b`, etc.

# Examples
```jldoctest
julia> u, v = Multivector.(basis(2))
2-element Vector{Multivector{2, 1, Vector{Int64}}}:
 1v1
 1v2

julia> using MacroTools: prettify

julia> ex = GeometricAlgebra.generated_multivector_function(*, u, v) |> prettify
:(let a = (CompositeMultivector(a)).components, b = (CompositeMultivector(b)).components
      comps = create_array(Vector{Any}, Int64, Val{1}(), Val{(4,)}(), a[1] * b[1] + a[2] * b[2], 0, 0, a[1] * b[2] + (-1 * a[2]) * b[1])
      (MixedMultivector{2})(comps)
  end)

```
"""
function generated_multivector_function(f, x...)
	syms = Symbol.('a' .+ (1:length(x)) .- 1)
	
	x_symb = symbolic_multivector.(x, syms)
	y_symb = f(x_symb...)
	comps = y_symb.components

	T = numberorany(promote_type(eltype.(x)...))
	comps_expr = SymbolicUtils.Code.MakeArray(comps, typeof(comps), T)

	assignments = [:( $sym = CompositeMultivector($sym).components ) for sym in syms]
	quote
		let $(assignments...)
			comps = $(SymbolicUtils.Code.toexpr(comps_expr))
			$(constructor(y_symb))(comps)
		end
	end
end
