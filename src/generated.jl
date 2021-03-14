index2ublade(::Type{<:Multivector{sig,k}}, i) where {sig,k} = lindex2ublade(UInt, k, i)
index2ublade(::Type{<:MixedMultivector}, i) = unsigned(i - 1)
ncomponents(::Type{<:Multivector{sig,k}}) where {sig,k} = binomial(dimension(sig), k) 
ncomponents(::Type{<:MixedMultivector{sig}}) where {sig} = 2^dimension(sig)

termexpr(s, i, j) = isone(s) ? :( a[$i]*b[$j] ) : isone(-s) ? :( -(a[$i]*b[$j]) ) : :( $s*a[$i]*b[$j] )

CM_AV = CompositeMultivector{<:AbstractVector}
function geometric_prod_gen_impl(A::Type{<:CM_AV}, B::Type{<:CM_AV})
	T = best_type(MixedMultivector, A, B)
	sig = signature(T)

	# quoted vector of `+(a, b, ...)` expressions
	iden = zero(eltype(T))
	comps = :([ $([:( +($iden) ) for _ ∈ 1:2^dimension(T)]...) ])

	# precompute multiplication table
	for i ∈ 1:ncomponents(A), j ∈ 1:ncomponents(B)
		s, u = ubladeprod(sig, index2ublade(A, i), index2ublade(B, j))
		iszero(s) && continue
		terms = comps.args[begin + u].args
		push!(terms, termexpr(s, i, j))
	end

	quote
		a, b = A.comps, B.comps
		$T($comps)
	end
end

function homogeneous_prod_gen_impl(A::Type{<:CM_AV}, B::Type{<:CM_AV}, k)
	T = best_type(Multivector, A, B; grade=Val(k))
	sig = signature(T)

	# quoted vector of `+(a, b, ...)` expressions
	iden = zero(eltype(T))
	comps = :([ $([:( +($iden) ) for _ ∈ 1:binomial(dimension(T), k)]...) ])

	# precompute multiplication table
	for i ∈ 1:ncomponents(A), j ∈ 1:ncomponents(B)
		s, u = ubladeprod(sig, index2ublade(A, i), index2ublade(B, j))
		if iszero(s) || ublade_grade(u) != k
			continue
		end
		terms = comps.args[ublade2lindex(u)].args
		push!(terms, termexpr(s, i, j))
	end

	quote
		a, b = A.comps, B.comps
		$T($comps)
	end
end

@generated geometric_prod_gen(A::CM_AV, B::CM_AV) = geometric_prod_gen_impl(A, B)
@generated homogeneous_prod_gen(A::CM_AV, B::CM_AV, ::Val{k}) where k = homogeneous_prod_gen_impl(A, B, k)