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

parse_sig(sig::AbstractString) = Tuple(get(Dict('+'=>1, '-'=>-1, '0'=>0), i) do i
	error("Invalid character in signature string: $i")
end for i ∈ sig)

function parse_sig(ctx, args...)
	if length(args) == 1
		# signature specified by literal
		sig = @eval ctx $(first(args))
		if sig isa AbstractString
			return parse_sig(sig)
		end
		return sig
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
			error("Invalid syntax $bv in signature")
		end
		push!(labels, label)
		push!(squares, sq)
	end
	NamedTuple{tuple(labels...)}(squares)
end

macro basis(args...)
	generate_blades(powerset, parse_sig(__module__, args...))
end

macro basisfull(args...)
	generate_blades(parse_sig(__module__, args...)) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end