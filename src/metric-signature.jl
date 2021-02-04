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

	dimension(::Sig) # if known. If omitted, treated as arbitrary/infinite-dimensional
	signature_labels(::Sig)
	show_signature(::Sig) # for pretty printing
=#


#= BASIS VECTOR LABELS =#

DEFAULT_BASIS_SYMBOL = "v"

const subscript_nums = '₀':'₉'
subscriptnum(n::Integer) = join(subscript_nums[begin + i] for i ∈ reverse(digits(n)))

signature_label(sig, i) = Symbol("$DEFAULT_BASIS_SYMBOL$(subscriptnum(i))")
signature_label(sig::NamedTuple, i::Integer) = keys(sig)[i]
signature_label(sig, i::Symbol) = i
signature_label(::MetricSignature{sig}, i) where sig = signature_label(sig, i)
signature_label(::OffsetSignature{sig}, i) where sig = signature_label(sig, i)

signature_labels(sig) = [signature_label(sig, i) for i ∈ 1:dimension(sig)]


# defaults & fallbacks

# signature_labels(sig) = error("cannot access components of metric signature $sig by label")
show_signature(sig) = repr(sig) # fallback
show_signature(sig::Type) = nameof(sig)

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

Base.getindex(::OffsetSignature{sig}, i) where sig = sig[i]
dimension(::OffsetSignature{sig}) where sig = dimension(sig)
show_signature(::OffsetSignature{sig,indices}) where {sig,indices} = "$(show_signature(sig))[$indices]"



unoffset_index(sig, i) = i
unoffset_index(sig, i::Symbol) = findfirst(==(i), signature_labels(sig))
unoffset_index(sig::Union{Tuple,NamedTuple}, i) = let r = 1:dimension(sig)
	i ∈ r ? i : error("index $i is outside range $r for signature $sig")
end
function unoffset_index(osig::OffsetSignature{sig,indices}, ioffset) where {sig,indices}
	i = findfirst(==(ioffset), indices)
	isnothing(i) && error("index $ioffset is outside range $indices for OffsetSignature $osig")
	i
end


sig_has_dim(sig) = applicable(dimension, sig)
sig_has_dim(::Type{EuclideanSignature}) = false


Minkowski = OffsetSignature{(t=-1,x=1,y=1,z=1),0:3}()



#= PRETTY PRINTING =#

Base.show(io::IO, sig::AbstractMetricSignature) = print(io, show_signature(sig))
function Base.show(io::IO, ::MIME"text/plain", sig::AbstractMetricSignature)
	print(io, "$sig = ")
	Base.show_default(io, sig)
end

# overload how AbstractMultivector types are displayed, in order to pretty-print metric signatures
function Base.show(io::IO, T::Type{<:AbstractMultivector})
	if isconcretetype(T)
		name = nameof(T)
		pretty_sig = show_signature(signature(T))
		params = T.parameters[2:end]
		print(io, "$name{$pretty_sig, $(join(params, ", "))}")
	else
		invoke(show, Tuple{IO,Type}, io, T) # call original show method
	end
end

