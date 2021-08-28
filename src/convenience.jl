basis_blades(sig, grade=1) = [Blade{sig}(1, i) for i ∈ bits_of_grade(grade, dimension(sig))]

basis(sig; grade=1) = basis_blades(sig, grade)
basis(dim::Integer; grade=1) = basis(Tuple(fill(1, dim)); grade=1)
basis(osig::OffsetSignature{sig,indices}; grade=1) where {sig,indices} = OffsetArray(basis_blades(osig, grade), indices)




# return expression containing variable assignments for each basis blade for the given signature.
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

parse_sig(n::Integer) = Tuple(ones(n))
parse_sig(sig::MetricSignature) = sig
parse_sig(sig::Tuple{Vararg{<:Integer}}) = sig
parse_sig(sig::NamedTuple) = sig
parse_sig(labels::Tuple{Vararg{Symbol}}) = NamedTuple{labels}(ones(Int, length(labels)))
parse_sig(sig::AbstractString) = Tuple(get(Dict('+'=>1, '-'=>-1, '0'=>0), i) do i
	throw(ArgumentError("Invalid character in signature string: $i"))
end for i ∈ sig)
parse_sig(sig) = throw(ArgumentError("could not interpret $sig as metric signature"))

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
			throw(ArgumentError("Invalid term '$bv' in signature specification"))
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
the label and square (defaulting to `1`) of a basis vector.
For an `n`-dimensional signature, `2^n` basis blades are defined.

See also `@basisall`.

Examples
===
```jldoctest
julia> @basis "+++"
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> @basis (+1,0,-1)
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> @basis x y
[ Info: Defined basis blades x, y, xy

julia> @basis dt=-1 dx=1 dy
[ Info: Defined basis blades dt, dx, dy, dtdx, dtdy, dxdy, dtdxdy

julia> dt^2 + dxdy*dy
MixedMultivector{⟨dt-,dx+,dy+⟩, Vector{Int64}}:
 -1
 1 dx
```
"""
macro basis(a, args...)
	generate_blades(powerset, parse_sig(__module__, a, args...))
end

"""
	@basisall

Similarly to `@basis`, populate the local namespace with basis blades,
but including all permutations of the basis blades.

Note that more than `2^n` variables are defined for a signature of dimension `n`!

Example
===
```jldoctest
julia> @basisall x y z
[ Info: Defined basis blades x, y, z, xy, yx, xz, zx, yz, zy, xyz, xzy, yxz, yzx, zxy, zyx

julia> zyx
Grade-3 Blade{⟨x+,y+,z+⟩, 3, Int64}:
 -1 xyz
```
"""
macro basisall(a, args...)
	generate_blades(parse_sig(__module__, a, args...)) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end