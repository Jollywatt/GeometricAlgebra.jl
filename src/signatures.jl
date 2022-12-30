#= Metric Signature Interface

Geometric algebras are defined by a metric signature, which
is any `isbitstype` object implementing:

- `dimension` giving the dimension of the underlying vector space
- `basis_vector_norm` for getting the norm of the ``i``th orthonormal basis vector

The first type parameter of an `AbstractMultivector` subtype is the
algebra’s defining signature.

=#


"""
	ncomponents(sig)
	ncomponents(sig, k)

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
⟨+---⟩

julia> BasisBlade{sig}
BasisBlade{⟨+---⟩} (pretty-printed BasisBlade{(1, -1, -1, -1)})
```
"""
show_signature(io, sig) = show(io, sig)
show_signature(io, sig::Tuple) = print(io, "⟨$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))⟩")


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
function show_basis_blade(io, sig, indices)
	if dimension(sig) < 10
		printstyled(io, "v"*join(string.(indices)); bold=true)
	else
		printstyled(io, join(string.("v", indices)); bold=true)
	end
end
show_basis_blade(io, sig::NamedTuple, indices) = printstyled(io, join(keys(sig)[indices]), bold=true)



"""
	componentstype(sig, N, T)

Array type to use to store components of multivectors of signature `sig`.
The resulting type should be able to store `N` components (in the case
of a fixed-size array) of type `T`.

The fallback method returns `Vector{T}` for `dimension(sig) <= 8`, and
`SparseVector{T}` otherwise.
"""
componentstype(sig, N, T) = dimension(sig) <= 8 ? Vector{T} : SparseVector{T}


#= Built-in Metric Signatures =#

# Int: Euclidean space
dimension(dim::Integer) = dim
basis_vector_norm(::Integer, i) = 1

# Tuple, NamedTuple
dimension(sig::Union{Tuple,NamedTuple}) = length(sig)
basis_vector_norm(sig::Union{Tuple,NamedTuple}, i) = sig[i]

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
dimension(::Cl{P,Q,R}) where {P,Q,R} = P + Q + R
basis_vector_norm(::Cl{P,Q,R}, i) where {P,Q,R} = i <= P ? 1 : i <= P + Q ? -1 : 0
show_signature(io, ::Cl{P,Q,R}) where {P,Q,R} = print(io, "Cl($P,$Q,$R)")
show_signature(io, ::Cl{P,Q,0}) where {P,Q} = print(io, "Cl($P,$Q)")
Base.show(io::IO, ::MIME"text/plain", sig::Cl) = show_pretty(io, show_signature, sig)


# experimental
struct SigWithStorageType{Sig,S} end
dimension(::SigWithStorageType{Sig}) where {Sig} = dimension(Sig)
basis_vector_norm(::SigWithStorageType{Sig}, i) where {Sig} = basis_vector_norm(Sig, i)
componentstype(::SigWithStorageType{Sig,S}, N, T) where {Sig,S<:StaticVector} = T <: Number ? S{N,T} : Vector{T}
componentstype(::SigWithStorageType{Sig,<:MVector}, N, T) where {Sig} = isbitstype(T) ? MVector{N,T} : Vector{T} # MVectors only support setindex! for isbits types
componentstype(::SigWithStorageType{Sig,SparseVector}, N, T) where {Sig} = SparseVector{T}
function show_signature(io, ::SigWithStorageType{Sig,S}) where {Sig,S}
	show_signature(io, Sig)
	print(io, " with ", nameof(S))
end



#= Convenience =#

interpret_signature(sig::String) = Tuple(Dict('+' => +1, '-' => -1, '0' => 0)[i] for i in sig)
interpret_signature(sig) = sig

"""
	basis(sig; grade=1)

Vector of basis blades of specified grade(s) for the geometric algebra defined by the metric signature `sig`.
The value `grade=:all` is a shortcut for `grade=0:dimension(sig)`.

See also [`@basis`](@ref) and [`@basisall`](@ref).

# Examples
```jldoctest
julia> basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 v1
 v2
 v3

julia> basis("-+++", grade=0:2:4)
8-element Vector{BasisBlade{⟨-+++⟩, _A, Int64} where _A}:
 1
 v12
 v13
 v23
 v14
 v24
 v34
 v1234

julia> basis(Cl(1,3), grade=:all) |> sum
16-component Multivector{Cl(1,3), 0:4, Vector{Int64}}:
 1
 1 v1 + 1 v2 + 1 v3 + 1 v4
 1 v12 + 1 v13 + 1 v23 + 1 v14 + 1 v24 + 1 v34
 1 v123 + 1 v124 + 1 v134 + 1 v234
 1 v1234
```
"""
function basis(sig; grade=1)
	sig = interpret_signature(sig)
	dim = dimension(sig)
	grade == :all && (grade = 0:dimension(sig))
	bits = componentbits(Val(dim), Val(grade))
	BasisBlade{sig}.(bits .=> 1)
end

# make expr assigning each of combos(basis(sig)) to a variable
function generate_blades(combos, sig)
	bvs = basis(sig)
	defs = Pair{Symbol}[]
	for (ordered_bvs, indices) ∈ zip(combos(bvs), combos(1:dimension(sig)))
		varname = sprint(show_basis_blade, sig, indices)
		isempty(varname) && continue
		push!(defs, Symbol(varname) => prod(ordered_bvs))
	end
	quote
		@info "Defined basis blades $($(join(first.(defs), ", ")))"
		$([ :($(esc(k)) = $v) for (k, v) ∈ defs ]...)
		nothing
	end
end

"""
	@basis sig

Populate namespace with basis blades of every grade in the geometric
algebra with metric signature `sig`.

!!! warning
	This defines ``2^n`` variables for an ``n`` dimensional signature.

See also [`@basisall`](@ref).

# Examples

```jldoctest
julia> @basis 3
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> 1v2 + 3v12
8-component Multivector{3, 0:3, Vector{Int64}}:
 1 v2
 3 v12
```
"""
macro basis(sig)
	generate_blades(powerset, interpret_signature(Main.eval(sig)))
end

"""
	@basisall sig

Similarly to [`@basis`](@ref), populate namespace with basis blades, but
include all permutations of each blade (e.g., `v21` as well as `v12`).

!!! warning
	This defines more than ``2^n`` variables for an ``n`` dimensional signature!


# Examples

```jldoctest
julia> @basisall (+1,-1)
[ Info: Defined basis blades v, v1, v2, v12, v21

julia> v12 == -v21
true
```
"""
macro basisall(sig)
	generate_blades(interpret_signature(Main.eval(sig))) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end


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

julia> cayleytable(basis((t=-1, x=1, y=1, z=1), grade=2), ∧)
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
	cayleytable(io, basis(sig, grade=:all), args...; kwargs...)
end

function cayleytable(io::IO, mvs::AbstractVector, op=*; separators=:auto, title=:( $(nameof(op))($(Symbol("↓")), $(Symbol("→"))) ))
	table = Any[op(a, b) for a ∈ mvs, b ∈ mvs]
	mvs_str = string.(mvs)

	if separators == :auto
		types = typeof.(mvs)
		diffs = types[begin + 1:end] .!= types[begin:end - 1]
		separators = [1; 1 .+ findall(diffs)]
	end

	pretty_table(
		io,
		string.(table),
		header = mvs_str,
		row_labels = mvs_str,
		row_label_column_title = string(title),
		vlines = separators,
		hlines = separators,
		crop = :horizontal,
	)
end