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

julia> Blade{sig}
Blade{⟨+---⟩} (pretty-printed Blade{(1, -1, -1, -1)})
```
"""
show_signature(io, sig) = show(io, sig)
show_signature(io, sig::Tuple) = print(io, "⟨$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))⟩")


"""
	show_basis_blade(io, sig, indices::Vector{Int})

Show the basis blade with unit vectors in `indices` for the geometric algebra
defined by `sig`.
Methods dispatching on `sig` should be added to customise basis blade labels
for particular algebras.

# Examples
```julia
julia> GeometricAlgebra.show_basis_blade(stdout, (1, 1, 1), [1, 3])
v13

julia> using GeometricAlgebra: subscript

julia> GeometricAlgebra.show_basis_blade(io, sig, indices) = print(io, join("𝒆".*subscript.(indices), "∧"))

julia> prod(basis(4))
Blade{⟨++++⟩, 4, Int64} of grade 4:
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

Default array type used to store components of multivectors of signature `sig`.
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

# Cl(p, q, r)
"""
	Cl(p, q=0, r=0)

Metric signature where `p`, `q` and `r` is the number of
basis vectors of norm `+1`, `-1` and `0`, respectively.

# Examples
```jldoctest
julia> basis(Cl(1,3))
4-element Vector{Blade{Cl(1,3), 1, Int64}}:
 1v1
 1v2
 1v3
 1v4

julia> ans .^ 2
4-element Vector{Blade{Cl(1,3), 0, Int64}}:
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
struct MetricWithStorage{Sig,S} end
dimension(::MetricWithStorage{Sig}) where {Sig} = dimension(Sig)
basis_vector_norm(::MetricWithStorage{Sig}, i) where {Sig} = basis_vector_norm(Sig, i)
componentstype(::MetricWithStorage{Sig,S}, N, T) where {Sig,S<:StaticVector} = S{N,T}
componentstype(::MetricWithStorage{Sig,SparseVector}, N, T) where {Sig} = SparseVector{T}
function show_signature(io, ::MetricWithStorage{Sig,S}) where {Sig,S}
	show_signature(io, Sig)
	print(io, " with $(nameof(S))")
end



#= Convenience =#

interpret_signature(sig::String) = Tuple(Dict('+' => +1, '-' => -1, '0' => 0)[i] for i in sig)
interpret_signature(sig) = sig

"""
	basis(sig; grade=1)

Return basis blades for the geometric algebra defined by the metric signature `sig`.

See also [`@basis`](@ref) and [`@basisall`](@ref).

# Examples
```jldoctest
julia> basis(3)
3-element Vector{Blade{3, 1, Int64}}:
 1v1
 1v2
 1v3

julia> prod(basis("-+++"))

julia> basis(Cl(1,3); grade=2)
6-element Vector{Blade{Cl(1,3), 2, Int64}}:
 1v12
 1v13
 1v23
 1v14
 1v24
 1v34
```
"""
basis(sig; grade=1) = let sig = interpret_signature(sig)
	if grade == :all
		Blade{sig}.(mmv_bits(Val(dimension(sig))) .=> 1)
	else
		Blade{sig}.(bits_of_grade(grade, dimension(sig)) .=> 1)
	end
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

See also [`@basisall`](@ref).

# Examples

```jldoctest
julia> @basis 3
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> 1v2 + 3v12
8-component MixedMultivector{3, Vector{Int64}}:
 1 v2
 3 v12
```
"""
macro basis(sig)
	generate_blades(powerset, interpret_signature(eval(sig)))
end

"""
	@basisall sig

Similarly to [`@basis`](@ref), populate namespace with basis blades, but
include all permutations of each blade.

!!! warning
	This defines ``2^n`` variables for an ``n`` dimensional signature!

# Examples

```jldoctest
julia> @basisall (+1,-1)
[ Info: Defined basis blades v, v1, v2, v12, v21

julia> v12 == -v21
true
```
"""
macro basisall(sig)
	generate_blades(interpret_signature(eval(sig))) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end


"""
	cayleytable(sig, op=*)

Display the multiplication table for a geometric algeba

The first argument may be a metric signature of a vector of elements
of the geometric algebra, and `op` may be any binary operator.

# Examples
```jldoctest
julia> cayleytable((t=-1, x=1, y=1))
 (↓) * (→) │   1 │   t     x    y │  tx    ty   xy │ txy
───────────┼─────┼────────────────┼────────────────┼─────
         1 │   1 │   t     x    y │  tx    ty   xy │ txy
───────────┼─────┼────────────────┼────────────────┼─────
         t │   t │  -1    tx   ty │  -x    -y  txy │ -xy
         x │   x │ -tx     1   xy │  -t  -txy    y │ -ty
         y │   y │ -ty   -xy    1 │ txy    -t   -x │  tx
───────────┼─────┼────────────────┼────────────────┼─────
        tx │  tx │   x     t  txy │   1    xy   ty │   y
        ty │  ty │   y  -txy    t │ -xy     1  -tx │  -x
        xy │  xy │ txy    -y    x │ -ty    tx   -1 │  -t
───────────┼─────┼────────────────┼────────────────┼─────
       txy │ txy │ -xy   -ty   tx │   y    -x   -t │   1

julia> cayleytable(basis(3), ∧)
 (↓) ∧ (→) │   v1    v2   v3
───────────┼─────────────────
        v1 │    0   v12  v13
        v2 │ -v12     0  v23
        v3 │ -v13  -v23    0

```
"""
function cayleytable(sig, args...)
	dim = dimension(sig)
	grade_slices = mmv_slice.(Val(dim), Val.(0:dim))
	separators = first.(grade_slices) #∪ [0, dim + 1]
	cayleytable(basis(sig, grade=:all), args...; separators)
end

function cayleytable(mvs::AbstractVector, op=*; separators=[1])
	table = [op(a, b) for a ∈ mvs, b ∈ mvs]
	mvs_str = string.(mvs)
	table
	pretty_table(
		string.(table),
		header = mvs_str,
		row_names = mvs_str,
		row_name_column_title = string(:( $(nameof(op))($(Symbol("↓")), $(Symbol("→"))) )),
		vlines = separators,
		hlines = separators,
	)
end