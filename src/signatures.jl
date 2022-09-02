#=
The Metric Signature Interface

Signatures should be all-bits values (so they can be the first type parameter of AbstractMultivectors)
and should implement:
 - dimension(sig)
 - getindex(sig, i) to give the norm of the ith basis vector
and optionally:
 - default_storage_type(sig) preferred AbstractVector type for components (e.g., SVector, Vector or SparseVector)
 - basis_blade_label(sig, indices) for custom printing basis blades (e.g., "e12" or "e₁∧e₂" or "dxdy")
 - show_signature(sig) for pretty-printing in type parameters (e.g., "𝒢(3)" or "⟨+++⟩")

=#

# fallback with (named) tuples in mind
dimension(sig) = length(sig)

# `default_storagetype` should choose types appropriately by taking into account
# the algebra's dimension for optimal memory useage and performance
default_storagetype(sig, T) = dimension(sig) >= 8 ? SparseVector{UInt,T} : Vector{T}


"""
	basis_blade_label(sig, indices[, labels])

Human-readable string representation of the basis blade specified by `indices`.

Basis blade is specified by `indices`, a vector of integers.
"""
basis_blade_label(sig, indices, labels=1:dimension(sig)) = "v$(join(labels[indices]))"
basis_blade_label(sig::NamedTuple{labels}, indices, _) where labels = join(labels[indices])

abstract type MetricSignature end

"""
	OffsetSignature(sig, indices)

Metric signature `sig` with offset `indices`, enabling non-standard indexing for multivectors.

Example
===
```jldoctest
julia> lorentzian = OffsetSignature((-1,1,1,1), 0:3) # zero-based indexing
⟨-+++⟩ with indices 0:3
(pretty-printed OffsetSignature{(-1, 1, 1, 1), 0:3}())

julia> (1:4)'basis(lorentzian) # construct spacetime 4-vector
Grade-1 Multivector{⟨-+++⟩ with indices 0:3, 1, Vector{Int64}}:
 1 v0
 2 v1
 3 v2
 4 v3

julia> ans[0] # get first ("time") component
1
```
"""
struct OffsetSignature{sig,indices} <: MetricSignature end
OffsetSignature(sig,indices) = OffsetSignature{sig,indices}()
Base.getindex(::OffsetSignature{sig}, args...) where sig = getindex(sig, args...)
dimension(::OffsetSignature{sig}) where sig = dimension(sig)
basis_blade_label(::OffsetSignature{sig,labels}, indices) where {sig,labels} = basis_blade_label(sig, indices, labels)


"""
Convert human-readable basis vector index or symbol into 1-based integer index
"""
normalize_bv_index(sig, i) = i
normalize_bv_index(sig, i::Symbol) = findfirst(==(string(i)), [basis_blade_label(sig, [i]) for i ∈ 1:dimension(sig)])
function normalize_bv_index(sig::Union{Tuple,NamedTuple}, i::Integer)
	r = 1:dimension(sig)
	i ∈ r ? i : error("index $i is outside range $r for signature $sig")
end
function normalize_bv_index(osig::OffsetSignature{sig,indices}, ioffset::Integer) where {sig,indices}
	i = findfirst(==(ioffset), indices)
	isnothing(i) && error("index $ioffset is outside range $indices for $osig")
	i
end


"""
Pretty-print metric signature in short, non-parseable form such as ⟨-+++⟩
"""
show_signature(sig) = repr(sig) # fallback
show_signature(sig::Tuple) = "⟨$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))⟩"
function show_signature(sig::NamedTuple)
	items = map(zip(keys(sig), sig)) do (label, square)
		s = get(Dict(+1=>"+", -1=>"-"), square, square)
		"$label$s"
	end
	"⟨$(join(items, ","))⟩"
end
function show_signature(::OffsetSignature{sig,indices}) where {sig,indices}
	"$(show_signature(sig)) with indices $indices"
end