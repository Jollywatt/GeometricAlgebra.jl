function generate_blades(combos, sig)
	bvs = collect(basis(sig)) # needed because `OffsetArrays` don't work with, e.g., `powerset`
	vars = Symbol[]
	vals = Blade{sig}[]
	for (ordered_bvs, ublade) ∈ zip(combos(bvs), combos(1:dimension(sig)))
		varname = basis_blade_label(sig, ublade)
		isempty(varname) && continue
		push!(vars, Symbol(varname))
		push!(vals, prod(ordered_bvs))
	end
	quote
		@info "Defined basis blades $(join($vars, ", "))"
		$([ :($(esc(k)) = $(v)) for (k, v) ∈ zip(vars, vals) ]...)
		nothing
	end
end
generate_blades(sig) = generate_blades(x -> [[i] for i ∈ x], sig)

function parse_sig(args...)
	if length(args) == 1
		# signature specified by literal
		return @eval $(first(args))
	end

	# signature specified by series of exprs `label[=square]`
	labels = Symbol[]
	squares = []
	for bv ∈ args
		if bv isa Symbol
			label, sq = bv, 1
		elseif bv isa Expr && bv.head == :(=)
			label, sq = bv.args
		else
			error("invalid syntax $bv in signature")
		end
		push!(labels, label)
		push!(squares, sq)
	end
	NamedTuple{tuple(labels...)}(squares)
end

macro basis(args...)
	generate_blades(powerset, parse_sig(args...))
end

macro basisperm(args...)
	generate_blades(parse_sig(args...)) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end