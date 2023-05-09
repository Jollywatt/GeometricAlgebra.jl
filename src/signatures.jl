#= Metric Signature Interface

Geometric algebras are defined by a metric signature, which
is any `isbitstype` object implementing:

- `dimension` giving the dimension of the underlying vector space
- `basis_vector_square` for getting the scalar square of the ``i``th orthonormal basis vector

The first type parameter of an `<:AbstractMultivector{Sig}` is the
algebra’s defining metric signature.

In addition to the required methods above, metric signatures may implement

- `show_signature(io, sig)` for custom printing of the `Sig` type parameter
- `show_basis_blade(io, sig, indices)` for custom basis blade styles (e.g., "dx ∧ dy", "e₁₂")
=#


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
show_basis_blade(io::IO, sig::NamedTuple, indices::Vector) = printstyled(io, join(keys(sig)[indices]), bold=true)

show_basis_blade(io::IO, sig, bits::Unsigned) = show_basis_blade(io, sig, bits_to_indices(sig, bits))


"""
	componentstype(sig, N, T)

Array type to use to store components of multivectors of signature `sig`.
The resulting type should be able to store `N` components (in the case
of a fixed-size array) of type `T`.

The fallback method returns `MVector{N,T}` for `dimension(sig) <= 8`, and
`Vector{T}` otherwise.
"""
componentstype(sig, N, T) = dimension(sig) <= 8 ? MVector{N,T} : Vector{T}


#= Built-in Metric Signatures =#

# Int: Euclidean space
dimension(dim::Integer) = dim
basis_vector_square(::Integer, i) = 1

# Tuple, NamedTuple
dimension(sig::Union{Tuple,NamedTuple}) = length(sig)
basis_vector_square(sig::Union{Tuple,NamedTuple}, i) = sig[i]

"""
	Cl(p, q=0, r=0)

Metric signature where `p`, `q` and `r` are the number of
basis vectors of norm `+1`, `-1` and `0`, respectively.

See also [`@Cl_str`](@ref).

# Examples
```jldoctest
julia> basis(Cl(1,3))
4-element Vector{BasisBlade{Cl(1,3), 1, Int64}}:
 1 v1
 1 v2
 1 v3
 1 v4

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
dimension(::Cl{P,Q,R}) where {P,Q,R} = P + Q + R
basis_vector_square(::Cl{P,Q,R}, i) where {P,Q,R} = i <= P ? 1 : i <= P + Q ? -1 : 0
show_signature(io, ::Cl{P,Q,R}) where {P,Q,R} = print(io, "Cl($P,$Q,$R)")
show_signature(io, ::Cl{P,Q,0}) where {P,Q} = print(io, "Cl($P,$Q)")
Base.show(io::IO, ::MIME"text/plain", sig::Cl) = show_pretty(io, show_signature, sig)




#= Convenience =#

"""
	Cl(sig::String) -> Tuple

Shorthand for a tuple specifying a metric signature, e.g., `Cl("-+++") === (-1, +1, +1, +1)`.
String may contain `'+'`, `'-'` and `'0'`.

# Example
```jldoctest
julia> Cl("+++") # 3D Euclidean metric signature
(1, 1, 1)

julia> basis(ans)
3-element Vector{BasisBlade{Cl("+++"), 1, Int64}}:
 1 v1
 1 v2
 1 v3
```
"""
Cl(s::String) = interpret_signature(s)

# depreciated?
macro Cl_str(s)
	interpret_signature(s)
end


# show_signature(io, sig::Tuple) = print(io, "Cl\"$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))\"")
show_signature(io, sig::Tuple) = print(io, "Cl(\"$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))\")")

interpret_signature(sig::String) = Tuple(Dict('+' => +1, '-' => -1, '0' => 0)[i] for i in sig)
interpret_signature(sig) = sig


"""
	cayleytable(sig, op=*)
	cayleytable(objs, op=*)

Display a multivector multiplication table.

The first argument may be a metric signature or any vector of objects
which can be combined with the binary operator `op`.

The keyword argument `title` sets the contents of the top-left cell.

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

julia> cayleytable(basis((t=-1, x=1, y=1, z=1), 2), ∧)
 (↓) ∧ (→) │   tx     ty    xy    tz     xz    yz
───────────┼──────────────────────────────────────
        tx │    0      0     0     0      0  txyz
        ty │    0      0     0     0  -txyz     0
        xy │    0      0     0  txyz      0     0
        tz │    0      0  txyz     0      0     0
        xz │    0  -txyz     0     0      0     0
        yz │ txyz      0     0     0      0     0

```
"""
cayleytable(args...; kwargs...) = cayleytable(stdout, args...; kwargs...)
function cayleytable(io::IO, sig, args...; kwargs...)
	sig = interpret_signature(sig)
	cayleytable(io, basis(sig, :all), args...; kwargs...)
end

function cayleytable(io::IO, mvs::AbstractVector, op=*; separators=:auto, title=:( $(nameof(op))($(Symbol("↓")), $(Symbol("→"))) ))
	table = [op(a, b) for a ∈ mvs, b ∈ mvs]

	if separators == :auto
		types = typeof.(mvs)
		diffs = types[begin + 1:end] .!= types[begin:end - 1]
		separators = [1; 1 .+ findall(diffs)]
	end

	pretty_table(
		io,
		table,
		header = mvs,
		row_labels = mvs,
		row_label_column_title = string(title),
		vlines = separators,
		hlines = separators,
		crop = :horizontal,
	)
end