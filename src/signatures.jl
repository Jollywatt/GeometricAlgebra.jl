#= Metric Signature Interface

Geometric algebras are defined by a metric signature, which
is any `isbitstype` object implementing:

- `dimension` giving the dimension of the underlying vector space
- `getindex` for getting the norm of the ``i``th orthonormal basis vector

The first type parameter of an `AbstractMultivector` is the
algebra’s defining signature.

=#

dimension(sig::Union{Tuple,NamedTuple}) = length(sig)

"""
	show_signature(io, sig)

Pretty-print the metric signature `sig`.

This is used to display the metric signature type parameter
in `AbstractMultivector` subtypes to reduce visual noise.
Methods may optionally be added for user-defined metric signatures,
in a similar fashion to `Base.show`.

Examples
========

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

Examples
========
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



"""
	componentstype(sig, N, T)

Array type to use the components for multivectors of signature `sig`.
The resulting type should be able to store `N` components (in the case
of a fixed-size array) of element type `T`.

This is used when converting a `Blade` into a `CompositeMultivector` to
determine the type of the components array.
"""
componentstype(sig, N, T) = Vector{T}


struct EuclideanMetric
	dim::Int
end
dimension(sig::EuclideanMetric) = sig.dim
Base.getindex(::EuclideanMetric, i) = 1
show_signature(io, sig::EuclideanMetric) = printstyled(io, sig.dim, "D", bold=true)

struct MMetric{Sig} end
componentstype(::MMetric, N, T) = MVector{N,T}
dimension(::MMetric{Sig}) where {Sig} = dimension(Sig)
Base.getindex(::MMetric{Sig}, i) where {Sig} = Sig[i]




#= Convenience =#

basis(sig) = Blade{sig}.(bits_of_grade(1, dimension(sig)) .=> 1)
basis(dim::Integer) = basis(ntuple(i -> 1, dim))


function generate_blades(combos, sig)
	bvs = basis(sig)
	sig = signature(eltype(bvs))
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

Examples
========

```jldoctest
julia> @basis 3
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123

julia> 1v2 + 3v12
8-component MixedMultivector{⟨+++⟩, Vector{Int64}}:
 1 v2
 3 v12
```
"""
macro basis(sig)
	generate_blades(powerset, eval(sig))
end

"""
	@basisall sig

Similarly to [`@basis`](@ref), populate namespace with basis blades, but
include all permutations of each blade.

!!! warning
	This defines ``2^n`` variables for an ``n`` dimensional signature!

Examples
========

```jldoctest
julia> @basisall (+1,-1)
[ Info: Defined basis blades v, v1, v2, v12, v21

julia> v12 == -v21
true
```
"""
macro basisall(sig)
	generate_blades(eval(sig)) do bvs
		Iterators.flatten(permutations.(powerset(bvs)))
	end
end