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



function compile_multivector_function(f, x...)
	syms = Symbol.("x".*string.(1:length(x)))
	x_symb = symbolic_multivector.(x, syms)
	y_symb = f(x_symb...)

	comps = y_symb.components
end

