@generated function _generated_geometric_prod(a::CompositeMultivector{Sig}, b::CompositeMultivector{Sig}) where {Sig}
	ncomps = ncomponents(Sig)

	compexprs = [Expr[] for _ ∈ 1:ncomps]
	for (ai, abits) ∈ enumerate(bitsof(a)), (bi, bbits) ∈ enumerate(bitsof(b))
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		i = bits_to_mmv_index(bits, dimension(Sig))
		push!(compexprs[i], :( $factor * (ac[$ai] * bc[$bi]) ))
	end

	T = promote_type(eltype(a), eltype(b))
	S = componentstype(Sig, ncomps, T)
	if ismutabletype(S)
		assignments = [:( mmv.components[$i] = +($(terms...)) )
			for (i, terms) in enumerate(compexprs) if !isempty(terms)]
		quote
			ac, bc = a.components, b.components
			mmv = zero(MixedMultivector{$Sig,$S})
			$(assignments...)
			return mmv
		end
	else
		args = [isempty(terms) ? zero(T) : :( +($(terms...)) )
			for terms in compexprs]
		quote
			ac, bc = a.components, b.components
			comps = $S($(args...))
			MixedMultivector{$Sig,$S}(comps)
		end
	end
end

_generated_geometric_prod(a::CompositeMultivector{Sig}, b::Blade{Sig}) where {Sig} = _generated_geometric_prod(a, Multivector(b))
_generated_geometric_prod(a::Blade{Sig}, b::CompositeMultivector{Sig}) where {Sig} = _generated_geometric_prod(Multivector(a), b)