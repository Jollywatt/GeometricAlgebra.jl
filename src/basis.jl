struct BasisDisplayStyle
	blades::Dict{UInt,Tuple{Int,Vector{Int}}}
	blade_order
end


BASIS_DISPLAY_STYLES = IdDict(
	# 3 => BasisDisplayStyle(
	# 	Dict(0b101 => (-1, [3, 1])),
	# 	Dict(),
	# )
)
DEFAULT_BASIS_DISPLAY_STYLE = BasisDisplayStyle(Dict(), Dict())
get_basis_display_style(sig) = get(BASIS_DISPLAY_STYLES, sig, DEFAULT_BASIS_DISPLAY_STYLE)

function get_basis_blade(style::BasisDisplayStyle, bits::Unsigned)
	if bits ∈ keys(style.blades)
		style.blades[bits]
	else
		(1, bits_to_indices(bits))
	end
end




"""
	basis(sig, k=1)

Vector of basis blades of specified grade(s) `k` for the geometric algebra defined by the metric signature `sig`.
The value `k=:all` is a shortcut for `0:dimension(sig)`.

See also [`@basis`](@ref).

# Examples
```jldoctest
julia> basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 1 v1
 1 v2
 1 v3

julia> basis("-+++", 0:2:4)
8-element Vector{BasisBlade{Cl"-+++", _A, Int64} where _A}:
 1
 1 v12
 1 v13
 1 v23
 1 v14
 1 v24
 1 v34
 1 v1234

julia> basis(Cl(1,3), :all) |> sum
16-component Multivector{Cl(1,3), 0:4, MVector{16, Int64}}:
 1
 1 v1 + 1 v2 + 1 v3 + 1 v4
 1 v12 + 1 v13 + 1 v23 + 1 v14 + 1 v24 + 1 v34
 1 v123 + 1 v124 + 1 v134 + 1 v234
 1 v1234
```
"""
function basis(sig, k=1)
	sig = interpret_signature(sig)
	dim = dimension(sig)
	k == :all && (k = 0:dim)
	bits = componentbits(dim, k)
	BasisBlade{sig}.(1, bits)
end

"""
	basis(::Type{<:Multivector})

Create a generator which iterates over basis elements of the given multivector type.

# Examples
```jldoctest
julia> basis(Multivector{3,1}) |> collect
3-element Vector{Multivector{3, 1, MVector{3, Int64}}}:
 v1
 v2
 v3

julia> basis(Multivector{2,0:2,Vector{Bool}}) |> collect
4-element Vector{Multivector{2, 0:2, Vector{Bool}}}:
 (1)
 (v1)
 (v2)
 (v12)

```
"""
function basis(M::Type{<:Multivector})
	n = ncomponents(M)
	(M(SingletonVector(true, i, n)) for i in 1:n)
end
basis(M::Type{Multivector{Sig,K}}) where {Sig,K} = basis(Multivector{Sig,K,componentstype(Sig, ncomponents(M), Int)})


function generate_blades(sig;
                         grades=:all,
                         allperms=false,
                         pseudoscalar=:I,
                         scalar=false,
                         prefix=nothing,
                         style=get_basis_display_style(sig))
	
	grades = grades == :all ? (0:dimension(sig)) : grades
	grades = scalar ? grades : grades ∩ (1:dimension(sig))

	bits = componentbits(dimension(sig), grades) 

	if allperms
		indices = map(permutations∘bits_to_indices, bits) |> Iterators.flatten |> collect
	else
		indices = last.(get_basis_blade.(Ref(style), bits))
	end

	basisvectors = basis(sig)
	labels = sprint.(show_basis_blade, Ref(sig), indices)
	if !isnothing(prefix)
		labels .= replace.(labels, r"^[^0-9]+"=>prefix)
	end
	labels = Symbol.(labels)

	basisblades = labels .=> prod.(getindex.(Ref(basisvectors), indices))
	filter!(basisblades) do (label, _)
		label != Symbol("")
	end

	if dimension(sig) ∈ grades && !isnothing(pseudoscalar)
		push!(basisblades, Symbol(pseudoscalar) => prod(basisvectors))
	end

	basisblades
end

"""
	@basis sig grades=:all scalar=false pseudoscalar=:I allperms=false prefix=nothing

Populate namespace with basis blades for the geometric
algebra defined by metric signature `sig`.

Variable names are generated with [`show_basis_blade()`](@ref).

# Keyword arguments

- `grades`: which grades to define basis blades for (default `:all`).
- `scalar`: whether to include the unit scalar blade (e.g., `v`).
- `pseudoscalar`: alias for unit pseudoscalar (default `:I`).
  `pseudoscalar=nothing` defines no alias.
- `allperms`: include all permutations of each basis blade (e.g., define `v21` as well as `v12`).
- `prefix`: prefix for basis blades names (`nothing` leaves default names unchanged).

!!! warning
	This defines `2^dimension(sig)` variables with `grades=:all`, and more with `allperms=true`!

# Examples

```jldoctest
julia> @basis 3
[ Info: Defined basis blades v1, v2, v3, v12, v13, v23, v123, I in Main

julia> 1v2 + 3v12
8-component Multivector{3, 0:3, MVector{8, Int64}}:
 1 v2
 3 v12

julia> @basis "0++" prefix=:e
[ Info: Defined basis blades e1, e2, e3, e12, e13, e23, e123, I in Main

julia> @basis 2 allperms=true scalar=true pseudoscalar=nothing
[ Info: Defined basis blades v, v1, v2, v12, v21 in Main

julia> @basis (t=1,x=-1,y=-1,z=-1) grades=2 allperms=true
[ Info: Defined basis blades tx, xt, ty, yt, xy, yx, tz, zt, xz, zx, yz, zy in Main
```
"""
macro basis(sig, args...)
	sig = interpret_signature(__module__.eval(sig))
	pairs = @eval generate_blades($sig; $(args...))
	assignments = map(pairs) do (label, blade)
		:($(esc(label)) = $blade)
	end
	message = isempty(pairs) ?
		"No basis blades defined!" :
		"Defined basis blades $(join(string.(first.(pairs)), ", ")) in $__module__"

	# detect evaluation in top-level
	# see https://discourse.julialang.org/t/is-there-a-way-to-determine-whether-code-is-toplevel/39939/3
	canary = gensym("canary")
	quote
		$(assignments...)
		$(esc(canary)) = true
		if Base.isdefined($__module__, $(QuoteNode(canary)))
			@info $message
		end
	end
end
