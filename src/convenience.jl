
basisvarname(sig::Tuple, ublade) = "$DEFAULT_BASIS_SYMBOL"*join(map(string, ublade))
basisvarname(sig::NamedTuple, ublade) = join([signature_labels(sig)[i] for i ∈ ublade])

macro basis(sig, options=:all)
	sig = eval(sig)
	@assert dimension(sig) < 10
	bvs = basis(sig)
	vars = Symbol[]
	vals = []
	combos = Dict(
		:vector => x -> [[i] for i ∈ x],
		:all => powerset,
		:allcombos => x -> Iterators.flatten(permutations.(powerset(x)))
	)[options]
	for (bladebvs, ublade) ∈ zip(combos(bvs), combos(1:dimension(sig)))
		varname = basisvarname(sig, ublade)
		isempty(varname) && continue
		push!(vars, Symbol(varname))
		push!(vals, prod(bladebvs))
	end
	@info "Defined basis blades $(join(vars, ", "))"
	quote
		$([ :($(esc(k)) = $(v)) for (k, v) ∈ zip(vars, vals) ]...)
		nothing
	end
end

#=
perhaps:

@basis -> vector
@basisfull -> all
@basisfullperm -> allcombos

=#