#= METRIC SIGNATURES

The signature parameter `sig` in `AbstractMultivector{sig}` can simply be a
(named) tuple, which is good for most needs. But I want to keep things generic
to keep the possibility of infinite dimensions open, so `sig` could be
anything which implements `getindex`.

E.g., in high dimensions, `sig = EuclideanSignature(n)` is less cumbersome than
an n-tuple of 1s. Furthermore, if the dimension is left unspecified with
`sig = EuclideanSignature`, then multivectors can be implemented with ublades
of type `Vector{Symbol}`, and there is no mathematical restriction on
the number of dimensions.

This also allows for a declare-when-needed mode of interaction:
julia> x = basis(EuclideanSignature, :x)
1-Blade{EuclideanSignature, Float64, Vector{Symbol}, 1}
 1.0 x

julia> x + basis(EuclideanSignature, :y)
1-Multivector{EuclideanSignature, Dict{Vector{Symbol}, Float64}, 1}
 1.0 y
 1.0 x

=#

#= METRIC SIGNATURE INTERFACE

Required methods:

	getindex(::Sig, ::Integer)

Optional methods:

	dimension(::Sig) # only if known. If omitted, treated as arbitrary/infinite-dimensional

Methods for pretty printing:

	show_signature(::Sig)
	basis_vector_label(::Sig, i)
	basis_blade_label(::Sig, i)
=#


# TUPLE METRIC SIGNATURES
# signature_labels(sig::NamedTuple) = keys(sig)

dimension(sig::Union{Tuple,NamedTuple}) = length(sig)

show_signature(sig::Tuple) = "⟨$(join(map(s -> get(Dict(1 => "+", -1 => "-"), s, s), sig)))⟩"
function show_signature(sig::NamedTuple{labels,<:Tuple{Vararg{<:Integer}}}) where labels
	items = map(zip(keys(sig), sig)) do (label, square)
		s = get(Dict(1 => "+", -1 => "-"), square, square)
		"$label$s"
	end
	"⟨$(join(items, ","))⟩"
end

show_signature(sig::Type) = repr(sig)



abstract type AbstractMetricSignature end


struct MetricSignature{sig} <: AbstractMetricSignature end
MetricSignature(sig) = MetricSignature{sig}()

Base.getindex(::MetricSignature{sig}, i) where sig = sig[i]
dimension(::MetricSignature{sig}) where sig = dimension(sig)
show_signature(::MetricSignature{sig}) where sig = show_signature(sig)


struct EuclideanSignature <: AbstractMetricSignature
	n::Int
end

Base.getindex(::EuclideanSignature, i) = 1
dimension(sig::EuclideanSignature) = sig.n
show_signature(sig::EuclideanSignature) = "⟨$(sig.n)+⟩"

# used as open signature of unspecified dimension
dimension(sig::Type{EuclideanSignature}) = error("open signature $sig does not have specified dimension")
Base.getindex(::Type{EuclideanSignature}, i) = 1


struct OffsetSignature{sig,indices} <: AbstractMetricSignature end
OffsetSignature(sig, indices) = OffsetSignature{sig,indices}()

Base.getindex(::OffsetSignature{sig}, i) where sig = sig[i]
dimension(::OffsetSignature{sig}) where sig = dimension(sig)
show_signature(::OffsetSignature{sig,indices}) where {sig,indices} = "$(show_signature(sig))[$indices]"
Base.keys(::OffsetSignature{sig,indices}) where {sig,indices} = indices

offset_index(::Any, i) = i
offset_index(::OffsetSignature{sig,indices}, i) where {sig,indices} = indices[i]

unoffset_index(sig, i) = i
unoffset_index(sig, i::Symbol) = findfirst(==(i), basis_vector_labels(sig))
unoffset_index(sig::Union{Tuple,NamedTuple}, i::Integer) = let r = 1:dimension(sig)
	i ∈ r ? i : error("index $i is outside range $r for signature $sig")
end
function unoffset_index(osig::OffsetSignature{sig,indices}, ioffset::Integer) where {sig,indices}
	i = findfirst(==(ioffset), indices)
	isnothing(i) && error("index $ioffset is outside range $indices for OffsetSignature $osig")
	i
end


sig_has_dim(sig) = applicable(dimension, sig)
sig_has_dim(::Type{EuclideanSignature}) = false


Minkowski = OffsetSignature((t=-1,x=1,y=1,z=1), 0:3)



#= BASIS VECTOR LABELS

Add a `basis_vector_label(sig, i)` method for custom labelling for custom signature types.
Add a `basis_blade_label(sig, bvs::Vector)` method for more control; defaults to concatenating 
basis vector labels.

TODO:
redesign labelling interface
Need one function

=#

basis_vector_label(sig::NamedTuple, i::Integer) = keys(sig)[i]
basis_vector_label(sig, i::Symbol) = i

basis_vector_labels(sig) = [basis_vector_label(sig, i) for i ∈ 1:dimension(sig)]


function basis_blade_label(sig::T, bvs) where T
	if hasmethod(basis_vector_label, Tuple{T,eltype(bvs)})
		join(basis_vector_label(sig, i) for i ∈ bvs)
	else
		"v$(join(bvs))"
	end
end

function basis_blade_label(M::OffsetSignature{sig}, bvs) where sig
	if hasmethod(basis_vector_label, Tuple{typeof(sig),eltype(bvs)})
		join(basis_vector_label(sig, i) for i ∈ bvs)
	else
		"v$(join(offset_index(M, i) for i ∈ bvs))"
	end
end


struct Quaternions <: GeometricAlgebra.AbstractMetricSignature end

Base.getindex(::Quaternions, i) = 1
dimension(::Quaternions) = 3
function basis_blade_label(::Quaternions, bvs)
	get(Dict(
		[2,3] => :i,
		[1,3] => :j,
		[1,2] => :k,
	), bvs, basis_blade_label((1,1,1), bvs))
end