#= Metric Signature Interface

Geometric algebras are defined by a metric signature, which
is any `isbitstype` implementing the method

- `canonical_signature(sig) -> Tuple`

The first type parameter of an `<:AbstractMultivector{Sig}` is the
algebra’s defining metric signature.

In addition to the required methods above, metric signatures may implement

- `dimension(sig) == length(canonical_signature(sig))`
- `basis_vector_square(sig, i) == canonical_signature(sig)[i]`
- `componentstype(sig, N[, T])` for signature or size dependent default component storage types
- `show_signature(io, sig)` for custom printing of the `Sig` type parameter
- `show_basis_blade(io, sig, indices)` for custom basis blade styles (e.g., "dx ∧ dy", "e₁₂")
=#

"""
	canonical_signature(sig) -> Tuple

Canonical tuple representation of a metric signature.

# Examples
```jldoctest
julia> Cl(1,3)
Cl(1,3) (pretty-printed Cl{1, 3, 0}())

julia> GeometricAlgebra.canonical_signature(ans)
(1, -1, -1, -1)
```
"""
function canonical_signature end

dimension(sig) = length(canonical_signature(sig))
basis_vector_square(sig, i::Integer) = canonical_signature(sig)[i]


"""
	componentstype(sig, N) -> Type{<:AbstractVector}

The component array type for `N`-component multivectors with signature `sig`.

You can redefine this method to customise the default array type.
The fallback method returns `MVector{N}` for `N <= 16`, and `Vector` otherwise.
"""
componentstype(sig, N) = N <= 32 ? SVector{N} : Vector
componentstype(sig, N, T) = typeintersect(componentstype(sig, N), AbstractVector{T})

componentstype(T::Type{<:Multivector}) = componentstype(signature(T), ncomponents(T))
componentstype(::Multivector{Sig,K,T}) where {Sig,K,T} = T


"""
	show_signature(io, sig)

Pretty-print the metric signature `sig`.

This is used to display the metric signature type parameter
in `AbstractMultivector` subtypes to reduce visual noise.
Methods may optionally be added for user-defined metric signatures,
in a similar fashion to `Base.show`.

# Examples

```jldoctest
julia> sig = (+1,-1,-1,-1)
(1, -1, -1, -1)

julia> GeometricAlgebra.show_signature(stdout, sig)
Cl("+---")

julia> BasisBlade{sig}
BasisBlade{Cl("+---")} (pretty-printed BasisBlade{(1, -1, -1, -1)})
```
"""
show_signature(io, sig) = show(io, sig)

show_basis_blade(io::IO, sig, bits::Unsigned) = show_basis_blade(io, sig, bits_to_indices(sig, bits))

"""
	show_basis_blade(io, sig, indices::Vector{Int})

Show the basis blade ``𝒗_{i₁}⋯𝒗_{iₖ}`` with each ``iⱼ`` in `indices` in the geometric algebra defined by `sig`.
Methods dispatching on `sig` should be added to customise basis blade labels
for particular algebras.

# Examples
```julia
julia> GeometricAlgebra.show_basis_blade(stdout, (1, 1, 1), [1, 3])
v13

julia> using GeometricAlgebra: subscript

julia> GeometricAlgebra.show_basis_blade(io, sig, indices) = print(io, join("𝒆".*subscript.(indices), "∧"))

julia> prod(basis(4))
BasisBlade{⟨++++⟩, 4, Int64} of grade 4:
 1 𝒆₁∧𝒆₂∧𝒆₃∧𝒆₄
```
"""
function show_basis_blade(io::IO, sig, indices::Vector)
	if dimension(sig) < 10
		printstyled(io, "v"*join(string.(indices)); bold=true)
	else
		printstyled(io, join(string.("v", indices)); bold=true)
	end
end



"""
	ncomponents(sig) = 2^dimension(sig)
	ncomponents(sig, k) = binomial(dimension(sig), k)

Dimension of (the grade-`k` subspace of) the geometric algebra of metric
signature `sig`, viewed as a vector space.

If the dimension of the _underlying_ vector space (see [`dimension`](@ref)) in ``n``, then the algebra
is ``2^n``-dimensional, and its grade-``k`` subspace ``\\binom{n}{k}``-dimensional.
"""
ncomponents(sig) = 1 << dimension(sig)  # << constant folds whereas 2^dim doesn't
ncomponents(sig, k) = binomial(dimension(sig), k)



#= Built-in Metric Signatures =#

# Int: Euclidean space
canonical_signature(dim::Integer) = ntuple(_ -> 1, dim)
dimension(dim::Integer) = dim
basis_vector_square(::Integer, i::Integer) = 1

# Tuple, NamedTuple
canonical_signature(sig::Union{Tuple,NamedTuple}) = Tuple(sig)
dimension(sig::Union{Tuple,NamedTuple}) = length(sig)
basis_vector_square(sig::Union{Tuple,NamedTuple}, i::Integer) = sig[i]
show_basis_blade(io::IO, sig::NamedTuple, indices::Vector) = printstyled(io, join(keys(sig)[indices]), bold=true)

"""
	Cl(p, q=0, r=0)

Metric signature where `p`, `q` and `r` are the number of
basis vectors of norm `+1`, `-1` and `0`, respectively.

# Examples
```jldoctest
julia> basis(Cl(1,3))
4-element Vector{BasisBlade{Cl(1,3), 1, Int64}}:
 v1
 v2
 v3
 v4

julia> ans .^ 2
4-element Vector{BasisBlade{Cl(1,3), 0, Int64}}:
  1
 -1
 -1
 -1
```
"""
struct Cl{P,Q,R} end
Cl(p::Integer, q::Integer=0, r::Integer=0) = Cl{p,q,r}()
canonical_signature(::Cl{P,Q,R}) where {P,Q,R} = ntuple(i -> i <= P ? 1 : i <= P + Q ? -1 : 0, P + Q + R)
dimension(::Cl{P,Q,R}) where {P,Q,R} = P + Q + R
basis_vector_square(::Cl{P,Q,R}, i::Integer) where {P,Q,R} = i <= P ? 1 : i <= P + Q ? -1 : 0
show_signature(io, ::Cl{P,Q,R}) where {P,Q,R} = print(io, "Cl($P,$Q,$R)")
show_signature(io, ::Cl{P,Q,0}) where {P,Q} = print(io, "Cl($P,$Q)")
Base.show(io::IO, ::MIME"text/plain", sig::Cl) = show_pretty(io, show_signature, sig)




#= Convenience =#

"""
	Cl(sig::String) -> Tuple

Shorthand for a tuple specifying a metric signature, e.g., `Cl("-+++") === (-1, +1, +1, +1)`.
String may contain `'+'`, `'-'` and `'0'`.

For readability, `AbstractMultivector` types with a tuple metric signature display the signature as `Cl("...")`.

# Examples
```jldoctest
julia> Cl("+++") # 3D Euclidean metric signature
(1, 1, 1)

julia> basis(ans)
3-element Vector{BasisBlade{Cl("+++"), 1, Int64}}:
 v1
 v2
 v3

julia> Multivector{(0,-1,1,1,1),2}
Multivector{Cl("0-+++"), 2} (pretty-printed Multivector{(0, -1, 1, 1, 1), 2})
```
"""
Cl(s::String) = interpret_signature(s)

show_signature(io, sig::Tuple) = print(io, "Cl(\"$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))\")")

interpret_signature(sig::String) = Tuple(Dict('+' => +1, '-' => -1, '0' => 0)[i] for i in sig)
interpret_signature(sig) = sig




"""
	cayleytable(objs::AbstractVector, op=*; title)

Display a Cayley table for the binary operation `op` on the objects in `objs`.

The keyword argument `title` sets the contents of the top-left cell.

# Examples
```jldoctest
julia> cayleytable([false, true], &)
 (↓) & (→) │ false   true
───────────┼──────────────
     false │ false  false
      true │ false   true

julia> cayleytable(0:5, (a, b) -> mod(a*b, 6), title="ab mod 6")
 ab mod 6 │ 0  1  2  3  4  5
──────────┼──────────────────
        0 │ 0  0  0  0  0  0
        1 │ 0  1  2  3  4  5
        2 │ 0  2  4  0  2  4
        3 │ 0  3  0  3  0  3
        4 │ 0  4  2  0  4  2
        5 │ 0  5  4  3  2  1
```
"""
function cayleytable(io::IO, objs::AbstractVector, op=*; separators=:auto, title=:( $(nameof(op))($(Symbol("↓")), $(Symbol("→"))) ))
	table = [op(a, b) for a ∈ objs, b ∈ objs]

	if separators == :auto
		types = typeof.(objs)
		diffs = types[begin + 1:end] .!= types[begin:end - 1]
		separators = findall(diffs)
	end

	PrettyTables.pretty_table(
		io,
		table,
		column_labels = objs,
		row_labels = objs,
		stubhead_label = string(title),
		table_format=PrettyTables.TextTableFormat(
	        horizontal_lines_at_data_rows=separators,
	        vertical_lines_at_data_columns=separators,
	        horizontal_line_at_beginning=false,
	        horizontal_line_after_data_rows=false,
	        vertical_line_at_beginning=false,
	        vertical_line_after_data_columns=false,
	    )
	)
end



"""
	cayleytable(sig, op=*)

Display a multivector multiplication table for the geometric algebra with signature `sig`.

# Examples
```jldoctest
julia> cayleytable(3)
 (↓) * (→) │    1 │   v1     v2    v3 │  v12    v13   v23 │ v123
───────────┼──────┼───────────────────┼───────────────────┼──────
         1 │    1 │   v1     v2    v3 │  v12    v13   v23 │ v123
───────────┼──────┼───────────────────┼───────────────────┼──────
        v1 │   v1 │    1    v12   v13 │   v2     v3  v123 │  v23
        v2 │   v2 │ -v12      1   v23 │  -v1  -v123    v3 │ -v13
        v3 │   v3 │ -v13   -v23     1 │ v123    -v1   -v2 │  v12
───────────┼──────┼───────────────────┼───────────────────┼──────
       v12 │  v12 │  -v2     v1  v123 │   -1   -v23   v13 │  -v3
       v13 │  v13 │  -v3  -v123    v1 │  v23     -1  -v12 │   v2
       v23 │  v23 │ v123    -v3    v2 │ -v13    v12    -1 │  -v1
───────────┼──────┼───────────────────┼───────────────────┼──────
      v123 │ v123 │  v23   -v13   v12 │  -v3     v2   -v1 │   -1

julia> cayleytable(Cl(2,0,1), lcontract)
 lcontract(↓, →) │ 1 │ v1  v2  v3 │ v12  v13  v23 │ v123
─────────────────┼───┼────────────┼───────────────┼──────
               1 │ 1 │ v1  v2  v3 │ v12  v13  v23 │ v123
─────────────────┼───┼────────────┼───────────────┼──────
              v1 │ 0 │  1   0   0 │  v2   v3    0 │  v23
              v2 │ 0 │  0   1   0 │ -v1    0   v3 │ -v13
              v3 │ 0 │  0   0   0 │   0    0    0 │    0
─────────────────┼───┼────────────┼───────────────┼──────
             v12 │ 0 │  0   0   0 │  -1    0    0 │  -v3
             v13 │ 0 │  0   0   0 │   0    0    0 │    0
             v23 │ 0 │  0   0   0 │   0    0    0 │    0
─────────────────┼───┼────────────┼───────────────┼──────
            v123 │ 0 │  0   0   0 │   0    0    0 │    0

```
"""
cayleytable(args...; kwargs...) = cayleytable(stdout, args...; kwargs...)
cayleytable(io::IO, sig, args...; kwargs...) = cayleytable(io, basis(interpret_signature(sig), :all), args...; kwargs...)
