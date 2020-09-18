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

This also allows for a kinda cute mode of interaction:
julia> x = basis(EuclideanSignature, :x)
1-Blade{EuclideanSignature, Float64, Vector{Symbol}, 1}
 1.0 x

julia> x + basis(EuclideanSignature, :y)
1-Multivector{EuclideanSignature, Dict{Vector{Symbol}, Float64}, 1}
 1.0 y
 1.0 x

=#



dim(sig::Union{Tuple,NamedTuple}) = length(sig)

show_signature(sig::Tuple) = "⟨$(join(map(signum, sig)))⟩"
show_signature(sig::NamedTuple) = show_signature(values(sig))
show_signature(sig::NamedTuple{labels,<:Tuple{Vararg{<:Integer}}}) where labels  = "⟨$(join(map(((l,s),) -> "$l$(signum(s))", zip(keys(sig), sig)), ","))⟩"

abstract type AbstractMetricSignature end
show_signature(sig::Type) = nameof(sig)


signum(i) = i == 1 ? '+' : i == -1 ? '-' : "$i"

struct MetricSignature{sig} <: AbstractMetricSignature end
MetricSignature(sig) = MetricSignature{sig}()
Base.getindex(::MetricSignature{sig}, i) where sig = sig[i]
dim(::MetricSignature{sig}) where sig = dim(sig)
show_signature(::MetricSignature{sig}) where sig = show_signature(sig)

struct EuclideanSignature <: AbstractMetricSignature
	n::Int
end
Base.getindex(::EuclideanSignature, i) = 1
dim(sig::EuclideanSignature) = sig.n
show_signature(sig::EuclideanSignature) = "⟨$(sig.n)+⟩"

# open signature of unspecified dimension
Base.getindex(::Type{EuclideanSignature}, i) = 1
dim(::Type{EuclideanSignature}) = error("open signatures do not have specified dimension")

struct OffsetSignature{sig,indices} <: AbstractMetricSignature end
Base.getindex(::OffsetSignature{sig}, i) where sig = sig[i]
dim(::OffsetSignature{sig}) where sig = dim(sig)
show_signature(::OffsetSignature{sig,indices}) where {sig,indices} = "$(show_signature(sig))[$indices]"

Minkowski = OffsetSignature{(t=-1,x=1,y=1,z=1),0:3}()