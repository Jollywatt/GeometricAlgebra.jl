dimension(sig::Union{Tuple,NamedTuple}) = length(sig)

"""
	componentstype(sig)

Default array type for the components vector of a multivector,
used when converting a `Blade` into a `CompositeMultivector`.
"""
componentstype(sig) = Vector
componentstype(sig, T, N) = with_size(with_eltype(componentstype(sig), T), N)