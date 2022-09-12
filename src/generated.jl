@generated function _generated_geometric_prod(a::CompositeMultivector{Sig}, b::CompositeMultivector{Sig}) where {Sig}
	ncomps = mmv_size(Sig)

	comps = [Expr[] for _ ∈ 1:ncomps]
	for (ai, abits) ∈ enumerate(bitsof(a)), (bi, bbits) ∈ enumerate(bitsof(b))
		factor, bits = geometric_prod_bits(Sig, abits, bbits)
		i = bits_to_mmv_index(bits, dimension(Sig))
		push!(comps[i], :( $factor * (ac[$ai] * bc[$bi]) ))
	end
	assignments = [:( mmv.components[$i] = +($(comps[i]...)) )
		for i in 1:ncomps if length(comps[i]) > 0]

	T = promote_type(eltype(a), eltype(b))
	S = componentstype(Sig, ncomps, T)
	quote
		ac, bc = a.components, b.components
		mmv = zero(MixedMultivector{$Sig,$S})
		$(assignments...)
		return mmv
	end
end

_generated_geometric_prod(a::CompositeMultivector{Sig}, b::Blade{Sig}) where {Sig} = _generated_geometric_prod(a, Multivector(b))
_generated_geometric_prod(a::Blade{Sig}, b::CompositeMultivector{Sig}) where {Sig} = _generated_geometric_prod(Multivector(a), b)