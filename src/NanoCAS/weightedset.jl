#=

A "weighted set" is a collection of unique items each with
an associated non-zero weight. An element is a member of
a weighted set if its associated weight is non-zero.

=#


struct WeightedSet{K,V} <: AbstractDict{K,V}
	elements::Dict{K,V}
	function WeightedSet{K,V}(ps::AbstractDict) where {K,V}
		new{K,V}(filter!(!iszero∘last, ps))
	end
end

Base.empty(::Type{WeightedSet{K,V}}) where {K,V} = WeightedSet{K,V}(Dict{K,V}())

# when converting from types that might contain duplicate keys, accumulate repeats
function WeightedSet{K,V}(iter::Union{Tuple,AbstractVector,Base.Generator}) where {K,V}
	d = Dict{K,V}()
	for (k, v) in iter
		k in keys(d) ? d[k] += v : d[k] = v
	end
	WeightedSet{K,V}(d)
end

WeightedSet(d::Union{AbstractVector{Pair{K,V}},AbstractDict{K,V}}) where {K,V} = WeightedSet{K,V}(d)
function WeightedSet(ps::Pair...)
	K = promote_type(typeof.(first.(ps))...)
	V = promote_type(typeof.(last.(ps))...)
	WeightedSet{K,V}(ps)
end


pairtypes(::Type{Pair{L,R}}) where {L,R} = (L, R)
function WeightedSet(iter::Base.Generator)
	K, V = pairtypes(Base.@default_eltype(iter))
	WeightedSet{K,V}(iter)
end

Base.hash(a::WeightedSet, seed::UInt64) = hash(a.elements, seed)
Base.:(==)(a::WeightedSet, b::WeightedSet) = a.elements == b.elements

Base.copy(a::WeightedSet) = WeightedSet(copy(a.elements))

Base.iterate(a::WeightedSet) = iterate(a.elements)
Base.iterate(a::WeightedSet, state) = iterate(a.elements, state)
Base.length(a::WeightedSet) = length(a.elements)

Base.getindex(a::WeightedSet, k) = get(a.elements, k, zero(valtype(a)))
Base.setindex!(a::WeightedSet, v, k) = iszero(v) ? delete!(a.elements, k) : setindex!(a.elements, v, k)

combine(a::WeightedSet, b::WeightedSet) = WeightedSet(mergewith(+, a.elements, b.elements))
combine!(a::WeightedSet, b::WeightedSet) = WeightedSet(filter!(!iszero∘last, mergewith!(+, a.elements, b.elements)))

function mul!(a::WeightedSet, λ::Number)
	for k in keys(a)
		a[k] *= λ
	end
	a
end
negate!(a::WeightedSet) = mul!(a, -1)

