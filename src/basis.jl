"""
	BasisDisplayStyle(dim, blades[, blade_order]; kwargs...)
	BasisDisplayStyle(dim, blades_and_order; kwargs...)

Specifies how basis blades are displayed and ordered.
The default style for multivectors of metric signature `sig` can be set with
`GeometricAlgebra.BASIS_DISPLAY_STYLES[sig] = style`.

- `dim::Int` is the dimension of the algebra (number of basis vectors).
- `blades::Dict{UInt,Vector{Int}}` encodes the order of basis vectors
   in basis blades. E.g., `0b101 => [1, 3]` is the default style.
- `blade_order::Dict{Int,Vector{UInt}}` specifies the order of basis blades
   in a single grade. E.g., `3 => [0b011, 0b101, 0b110]` is the default ordering.
- `blades_and_order::Dict{Int,Vector{Int}}` gives a way of specifying the previous
   two mappings at once. E.g., `3 => [[1,2], [1,3], [2,3]]`.

# Keyword arguments

- `indices=1:dim` specifies the symbols used for each basis vector.
- `prefix="v"` is the prefix string for basis blades (if `sep == nothing`) or for each
   basis vector.
- `sep=nothing` is a string (e.g., `"∧"`) to separate each basis vector in a blade.
   If `sep` is `nothing`, blades are shown as e.g., `v123`, whereas an empty string
   results in `v1v2v3`.
- `labels` is a dictionary allowing individual basis blades to be given custom labels.
   E.g., `[3,2] => "𝒊"` means `4v32` is displayed as `4𝒊` (so long as the order
   `0b110 => [3,2]` is also specified in the `blades` argument — otherwise it would
   display as the default `-4v23`).


!!! note
	`BasisDisplayStyle` only affects how multivectors are _displayed_.
	The actual internal layout of multivectors is never affected.
	However, the active style for `sig` can affect the value of `basis(sig)`.

# Examples

```jldoctest; setup = :(delete!(GeometricAlgebra.BASIS_DISPLAY_STYLES, Cl(0,3)))
julia> Multivector{Cl(0,3),2}([3, -2, 1])
3-component Multivector{Cl(0,3), 2, Vector{Int64}}:
  3 v12
 -2 v13
  1 v23

julia> cyclical_style = BasisDisplayStyle(
           3, Dict(2 => [[2,3], [3,1], [1,2]]);
           indices = "₁₂₃",
           prefix = "e",
           sep = "",
           labels = Dict([1,2,3] => "I"),
       );

julia> GeometricAlgebra.BASIS_DISPLAY_STYLES[Cl(0,3)] = cyclical_style;

julia> Multivector{Cl(0,3),2}([3, -2, 1])
3-component Multivector{Cl(0,3), 2, Vector{Int64}}:
 1 e₂e₃
 2 e₃e₁
 3 e₁e₂

julia> ans*rdual(ans) # pseudoscalar `e₁e₂e₃` displayed as `I`
4-component Multivector{Cl(0,3), 1:2:3, MVector{4, Int64}}:
 14 I
```
To recover the default style:
```jldoctest
julia> delete!(GeometricAlgebra.BASIS_DISPLAY_STYLES, Cl(0,3))
IdDict{Any, BasisDisplayStyle}() 
```
"""
struct BasisDisplayStyle
	dim::Int
	blades::Dict{UInt,Vector{Int}}
	parities::Dict{UInt,Bool}
	order::Dict{Int,Vector{UInt}}
	indices::Vector{String}
	prefix::String
	sep::Union{Nothing,String}
	labels::Dict{UInt,String}
	function BasisDisplayStyle(
			dim,
			blades,
			order=Dict{Int,Vector{UInt}}();
			indices=1:dim,
			prefix="v",
			sep=nothing,
			labels=Dict{UInt,String}()
		)

		for (bits, indices) ∈ blades
			@assert bits == indices_to_bits(indices) "$(sprint(niceshow, dim, bits)) == indices_to_bits($indices)"
		end

		for (k, bits) ∈ order
			@assert issetequal(componentbits(dim, k), bits) "Blade order for key $k must be a permutation of componentbits($dim, $k)"
		end

		parities = Dict(bits => parity(sortperm(indices)) for (bits, indices) ∈ blades)

		if keytype(labels) <: Vector{<:Integer}
			labels = Dict(indices_to_bits(indices) => label for (indices, label) ∈ labels)
		end

		indices = string.(collect(indices))
		@assert length(indices) == dim

		new(
			dim,
			blades,
			parities,
			order,
			indices,
			prefix,
			sep,
			labels,
		)
	end
end

function BasisDisplayStyle(dim, blades_and_order::Dict{<:Integer,<:Vector{<:Vector}}=Dict{Int,Vector{Vector{Int}}}(); kwargs...)
	all_indices = reduce(vcat, values(blades_and_order); init = Int[])
	blades = Dict(indices_to_bits.(all_indices) .=> all_indices)
	order = Dict(k => indices_to_bits.(indices) for (k, indices) ∈ blades_and_order)
	BasisDisplayStyle(dim, blades, order; kwargs...)
end

const BASIS_DISPLAY_STYLES = IdDict{Any,BasisDisplayStyle}()
get_basis_display_style(sig) = get(BASIS_DISPLAY_STYLES, sig, sig)


bits_to_indices(style, bits::Unsigned) = bits_to_indices(bits)
bits_to_indices(style::BasisDisplayStyle, bits::Unsigned) =
	bits ∈ keys(style.blades) ? style.blades[bits] : bits_to_indices(bits)


basis_blade_parity(style, bits) = false
basis_blade_parity(style::BasisDisplayStyle, bits::Unsigned) = get(style.parities, bits, false)


componentbits(style::BasisDisplayStyle, k::Integer) = k ∈ keys(style.order) ? style.order[k] : componentbits(style.dim, k)
componentbits(style::BasisDisplayStyle, K) = Iterators.flatten(Iterators.map(k -> componentbits(style, k), K))
componentbits(style::BasisDisplayStyle, a::OrType{Multivector}) = componentbits(style, grade(a))


function show_basis_blade(io::IO, style::BasisDisplayStyle, bits::Unsigned)
	bits ∈ keys(style.labels) && return printstyled(io, style.labels[bits], bold=true)
	show_basis_blade(io, style, bits_to_indices(style, bits))
end
function show_basis_blade(io::IO, style::BasisDisplayStyle, indices::Vector)
	str = IOBuffer()
	if isnothing(style.sep)
		print(str, style.prefix)
		join(str, style.indices[indices])
	else
		join(str, style.prefix.*style.indices[indices], style.sep)
	end
	printstyled(io, String(take!(str)); bold=true)
end


# show unsigned numbers as bitstrings as long as the algebras’ dimension
niceshow(io::IO, dim, n::Unsigned) = print(io, "0b", string(n, pad = dim, base = 2))
niceshow(io::IO, dim, n) = show(io, n)
niceshow(io::IO, dim, (a, b)::Pair) = (niceshow(io, dim, a); print(io, " => "); niceshow(io, dim, b))
function niceshow(io::IO, dim, a::Union{Vector,Dict})
	first = true
	print(io, "[")
	for i in a
		first || print(io, ", ")
		niceshow(io, dim, i)
		first = false
	end
	print(io, "]")
end


function Base.show(io::IO, style::BasisDisplayStyle)
	println(io, typeof(style), ":")

	for i ∈ fieldnames(BasisDisplayStyle)
		print(io, lpad(i, 9), ": ")
		field = getfield(style, i)
		niceshow(io, style.dim, field)
		println(io)
	end

	println(io, "  preview:")
	n = style.dim
	mat = fill("", n + 1, binomial(n, n÷2))
	for k in 0:n
		blades = basis(n, k; style)
		for (i, b) in enumerate(blades)
			mat[begin + k,i] = sprint(show_basis_blade, style, b.bits)
		end
	end

	pretty_table(io, [string.(0:n) mat];
		vlines=[1],
		hlines=[],
		show_header=false,
		columns_width=[7; fill(0, size(mat, 2))],
		alignment=[:r; fill(:l, size(mat, 2))],
	)
end




"""
	basis(sig, k=1)

Vector of basis blades of specified grade(s) `k` for the geometric algebra defined by the metric signature `sig`.
The value `k=:all` is a shortcut for `0:dimension(sig)`.

The particular basis blades returned by `basis` and their order reflects the signature’s `BasisDisplayStyle`.
You can guarantee the default style by using
```julia
	BasisBlade{sig}.(1, componentbits(dimension(sig), k))
```
instead of `basis(sig, k)`.

See also [`@basis`](@ref).

# Examples
```jldoctest
julia> basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 1 v1
 1 v2
 1 v3

julia> basis("-+++", 0:2:4)
8-element Vector{BasisBlade{Cl("-+++"), _A, Int64} where _A}:
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
function basis(sig, k=1; style=nothing)
	sig = interpret_signature(sig)
	dim = dimension(sig)
	k == :all && (k = 0:dim)
	isnothing(style) && (style = get_basis_display_style(sig))
	bits = componentbits(style, k)
	parities = basis_blade_parity.(Ref(style), bits)
	BasisBlade{sig}.((-1).^parities, bits)
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

	bits = componentbits(style, grades) 

	if allperms
		indices = map(permutations∘bits_to_indices, bits) |> Iterators.flatten |> collect
	else
		indices = bits_to_indices.(Ref(style), bits)
	end

	basisvectors = basis(sig)
	labels = sprint.(show_basis_blade, Ref(style), indices)
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
