function symbolic_components(label, dims...)
	var = SymbolicUtils.Sym{Array{length(dims),Real}}(label)
	indices = Iterators.product(Base.OneTo.(dims)...)
	[SymbolicUtils.Term{Real}(getindex, [var, I...]) for I in indices]
end



with_eltype(::Type{<:Multivector{Sig,K,C}}, T) where {Sig,K,C} = Multivector{Sig,K,with_eltype(C, T)}
with_eltype(::Type{<:MixedMultivector{Sig,C}}, T) where {Sig,C} = MixedMultivector{Sig,with_eltype(C, T)}

function symbolic_multivector(T::Type{<:CompositeMultivector{Sig,C}}, label) where {Sig,C}
	with_eltype(T, Any)(symbolic_components(label, ncomponents(T)))
end
symbolic_multivector(a::AbstractMultivector, label) = symbolic_multivector(typeof(a), label)


"""
	generated_multivector_function(f, a, b, ...)

Trace evaluation of `f(a, b, ...)` by evaluating `f` on symbolic versions of each
`CompositeMultivector` instance or type `a`, `b`, ..., returning an expression body
suitable for a `@generated` function.

The names of the arguments to be passed to the function
body are the literal symbols `:a`, `:b`, etc.

# Examples
```jldoctest
julia> u, v = Multivector.(basis(2))
2-element Vector{Multivector{⟨++⟩, 1, Vector{Int64}}}:
 1v1
 1v2

julia> using MacroTools: prettify

julia> ex = GeometricAlgebra.generated_multivector_function(*, u, v) |> prettify
:(let a = a.components, b = b.components
      (MixedMultivector{⟨++⟩})([a[1] * b[1] + a[2] * b[2], 0, 0, a[1] * b[2] + (-1 * a[2]) * b[1]])
  end)

```
"""
function generated_multivector_function(f, x...)
	syms = Symbol.('a' .+ (1:length(x)) .- 1)
	x_symb = symbolic_multivector.(x, syms)
	y_symb = f(x_symb...)

	comps = SymbolicUtils.Code.toexpr.(y_symb.components)
	assignments = [:( $sym = $sym.components ) for (i, sym) in enumerate(syms)]

	comps_expr = if comps isa Vector
		:( [$(comps...)] )
	elseif comps isa MVector
		:( MVector($(comps...)) )
	elseif comps isa SVector
		:( SVector($(comps...)) )
	end

	quote
		let $(assignments...)
			$(constructor(y_symb))($comps_expr)
		end
	end
end
