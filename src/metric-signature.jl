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



dim(sig::Union{Tuple,NamedTuple}) = length(sig)
signature_labels(sig::NamedTuple) = keys(sig)
signature_labels(sig) = error("cannot access components of metric signature $(show_signature(sig)) = $sig by label")

signum(i) = i == 1 ? '+' : i == -1 ? '-' : "$i"
show_signature(sig::Tuple) = "⟨$(join(map(signum, sig)))⟩"
show_signature(sig::NamedTuple) = show_signature(values(sig))
show_signature(sig::NamedTuple{labels,<:Tuple{Vararg{<:Integer}}}) where labels  = "⟨$(join(map(((l,s),) -> "$l$(signum(s))", zip(keys(sig), sig)), ","))⟩"
show_signature(sig::Type) = nameof(sig)



abstract type AbstractMetricSignature end

struct MetricSignature{sig} <: AbstractMetricSignature end
MetricSignature(sig) = MetricSignature{sig}()
dim(::MetricSignature{sig}) where sig = dim(sig)
Base.getindex(::MetricSignature{sig}, i) where sig = sig[i]
show_signature(::MetricSignature{sig}) where sig = show_signature(sig)


struct EuclideanSignature <: AbstractMetricSignature
	n::Int
end
dim(sig::EuclideanSignature) = sig.n
Base.getindex(::EuclideanSignature, i) = 1
show_signature(sig::EuclideanSignature) = "⟨$(sig.n)+⟩"

# used as open signature of unspecified dimension
dim(sig::Type{EuclideanSignature}) = error("open signature $sig does not have specified dimension")
Base.getindex(::Type{EuclideanSignature}, i) = 1


struct OffsetSignature{sig,indices} <: AbstractMetricSignature end
dim(::OffsetSignature{sig}) where sig = dim(sig)
Base.getindex(::OffsetSignature{sig}, i) where sig = sig[i]
show_signature(::OffsetSignature{sig,indices}) where {sig,indices} = "$(show_signature(sig))[$indices]"
signature_labels(::OffsetSignature{sig}) where sig = keys(sig)



unoffset_index(sig, i) = i
unoffset_index(sig, i::Symbol) = findfirst(==(i), signature_labels(sig))
unoffset_index(sig::Union{Tuple,NamedTuple}, i) = let r = 1:dim(sig)
	i ∈ r ? i : error("index $i is outside range $r for signature $sig")
end
function unoffset_index(osig::OffsetSignature{sig,indices}, ioffset) where {sig,indices}
	i = findfirst(==(ioffset), indices)
	isnothing(i) && error("index $ioffset is outside range $indices for OffsetSignature $osig")
	i
end


Minkowski = OffsetSignature{(t=-1,x=1,y=1,z=1),0:3}()

sig_has_dim(sig) = applicable(dim, sig)
sig_has_dim(::Type{EuclideanSignature}) = false