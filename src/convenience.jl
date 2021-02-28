"""
Return expression containing variable assignments for each basis blade for the given signature.
"""
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

parse_sig(sig::Tuple{Vararg{<:Integer}}) = sig
parse_sig(sig::NamedTuple) = sig
parse_sig(labels::Tuple{Vararg{Symbol}}) = NamedTuple{labels}(ones(Int, length(labels)))
parse_sig(sig::AbstractString) = Tuple(get(Dict('+'=>1, '-'=>-1, '0'=>0), i) do i
	error("Invalid character in signature string: $i")
end for i ∈ sig)
function parse_sig(sig::Tuple)
	# signature specified by series of symbols/exprs `label[=square]`
	labels = Symbol[]
	squares = []
	for bv ∈ sig
		if bv isa Symbol
			label, sq = bv, 1
		elseif bv isa Expr && bv.head == :(=)
			label, sq = bv.args
		else
			error("Invalid term '$bv' in signature specification")
		end
		push!(labels, label)
		push!(squares, sq)
	end
	NamedTuple{tuple(labels...)}(squares)
end

# for a single argument `@basis X`, interpret `X` as an expression which evaluates to a signature.
# for many arguments, `@basis x1 x2 x3...`, interpret `(x1, x2, x3)` as a signature.
parse_sig(::Module, args...) = parse_sig(args)
parse_sig(m::Module, args) = let a = @eval m $args
	parse_sig(a)
end

"""
	@basis sig
	@basis x y z...

Populate the local namespace with basis blades of all grades for
the geometric algebra with signature `sig`, which may be provided
as a string, tuple, named tuple or `<:AbstractMetricSignature`.
If multiple arguments are given, each is interpreted as specifying
a label and square (defaulting to `1`) of a basis vector.
For an `n`-dimensional signature, `2^n` basis blades are defined.

See also `@basis_perms`.

Examples
===
```jldoctest
julia> @basis "+++"
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> @basis (+1,0,-1)
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> @basis dt=-1 dx
[ Info: Defined basis blades dt, dx, dtdx

julia> dt^2 + dx
MixedMultivector{⟨dt-,dx+⟩, Array{Float64,1}}:
 -1.0
 1.0 dx
```
"""
macro basis(args...)
	generate_blades(powerset, parse_sig(__module__, args...))
end

"""
	@basis_perms

Similarly to `@basis`, populate the local namespace with basis blades,
but include names for all permutations of the basis vectors.

More than `2^n` variables are defined for a signature of dimension `n`!

Example
===
```jldoctest
julia> @basis_perms x y z
[ Info: Defined basis blades x, y, z, xy, yx, xz, zx, yz, zy, xyz, xzy, yxz, yzx, zxy, zyx

julia> zyx
Grade-3 Blade{⟨x+,y+,z+⟩, 3, Float64, UInt8}:
 -1.0 xyz
```
"""
macro basis_perms(args...)
	generate_blades(parse_sig(__module__, args...)) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end