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
	componentstype(sig)

Default array type for the components vector of a multivector,
used when converting a `Blade` into a `CompositeMultivector`.
"""
componentstype(sig) = Vector


struct EuclideanMetric
	dim::Int
end
dimension(sig::EuclideanMetric) = sig.dim
Base.getindex(::EuclideanMetric, i) = 1


#= Display Methods =#

basis_blade_label(sig, indices) = "v"*join(string.(indices))