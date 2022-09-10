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

struct SMetric
	dim::Int
end
componentstype(::SMetric, N, T) = MVector{N,T}
dimension(sig::SMetric) = sig.dim
Base.getindex(::SMetric, i) = 1

#= Display Methods =#

basis_blade_label(sig, indices) = "v"*join(string.(indices))