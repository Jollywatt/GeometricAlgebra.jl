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

If the dimension of the _underlying_ vector space in ``n``, then the algebra
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
	show_basis_blade(io, sig, indices)

Show the basis blade which is the product of the unit vectors in `indices`
in a geometric algebra defined by `sig`.
Methods dispatching on `sig` should be added to customise basis blade labels
for particular algebras.
	
The fallback method is:
```julia
show_basis_blade(io, sig, indices) = printstyled(io, "v"*join(string.(indices)); bold=true)
```

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
show_basis_blade(io, sig, indices) = printstyled(io, "v"*join(string.(indices)); bold=true)
show_basis_blade(io, sig::NamedTuple, indices) = printstyled(io, join(keys(sig)[indices]), bold=true)



"""
	componentstype(sig, N, T)

Array type to use the components for multivectors of signature `sig`.
The resulting type should be able to store `N` components (in the case
of a fixed-size array) of element type `T`.

This is used when converting a `Blade` into a `CompositeMultivector` to
determine the type of the components array.
"""
componentstype(sig, N, T) = Vector{T}


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
componentstype(::MetricWithStorage{Sig,S}, N, T) where {Sig,S} = S{N,T}
function show_signature(io, ::MetricWithStorage{Sig,S}) where {Sig,S}
	print(io, nameof(S))
	show_signature(io, Sig)
end



#= Convenience =#

interpret_signature(sig::String) = Tuple(Dict('+' => +1, '-' => -1, '0' => 0)[i] for i in sig)
interpret_signature(sig) = sig

"""
	basis(sig; grade=1)

Return basis blades for the geometric algebra defined by the metric signature `sig`.

# Examples
```jldoctest
julia> basis(3)
3-element Vector{Blade{3, 1, Int64}}:
 1v1
 1v2
 1v3

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
	Blade{sig}.(bits_of_grade(grade, dimension(sig)) .=> 1)
end

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