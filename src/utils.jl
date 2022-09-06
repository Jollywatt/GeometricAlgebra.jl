promote_to(T, x) = convert(promote_type(T, typeof(x)), x)

zeroslike(::Type{Vector{T}}, n) where {T} = zeros(T, n)

with_eltype(::Type{<:Vector}, T) = Vector{T}
with_size(V::Type{<:Vector}, N) = V