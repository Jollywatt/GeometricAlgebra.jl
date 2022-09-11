 @generated function _geometric_prod_gen(a::CompositeMultivector{Sig}, b::CompositeMultivector{Sig}) where {Sig}
	T = promote_type(eltype(a), eltype(b))
	comps = fill(:( +($(zero(T))) ), mmv_size(Sig))
	for (ai, abits) ∈ enumerate(bitsof(a)), (bi, bbits) ∈ enumerate(bitsof(b))
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		i = bits_to_mmv_index(bits, dimension(Sig))
		prev_terms = comps[i].args[2:end]
		comps[i] = :( +($(prev_terms...), $factor*(ac[$ai]*bc[$bi])) )
	end
	quote
		ac, bc = a.components, b.components
		comps = @inbounds $T[$(comps...)]
		MixedMultivector{$Sig}(comps)
	end
end

_geometric_prod_gen(a::CompositeMultivector{Sig}, b::Blade{Sig}) where {Sig} = _geometric_prod_gen(a, Multivector(b))
_geometric_prod_gen(a::Blade{Sig}, b::CompositeMultivector{Sig}) where {Sig} = _geometric_prod_gen(Multivector(a), b)