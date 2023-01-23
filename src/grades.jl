#= Grade promotion =#

"""
	promote_grades(dim, k)

Canonicalize the grade type parameter `k`.

Returns a subset of `0:dim`, while attempting to normalize equivalent
representations, such as `0:1:3 => 0:3` or `(3, 0) => (0, 3)`.
"""
promote_grades(dim::Integer, k::Integer) = k # its more helpful to have
function promote_grades(dim::Integer, k)
	k = (0:dim) ∩ k
	length(k) == 1 ? first(k) : Tuple(k)
end
function promote_grades(dim::Integer, k::OrdinalRange)
	lo, hi = max(0, minimum(k)), min(dim, maximum(k))
	lo == hi && return lo
	Δ = abs(step(k))
	Δ == 1 ? (lo:hi) : (lo:Δ:hi)
end

"""
	promote_grades(dim, p, q, ...)

Return a suitable grade type parameter which contains the grades `p ∪ q ∪ ...`.

In order to reduce the number of possible type parameters,
the result may be larger than the exact union.
Specifically, when combining different grades, `promote_grades` will try to return
the narrowest grade(s) out of:
 - an integer `k ∈ 0:dim` for homogeneous elements (fewest components)
 - `(0, dim)`, for elements in the scalar-pseudoscalar subalgebra
 - `0:2:dim`, for elements in the even subalgebra
 - `0:dim`, for general inhomogeneous elements (most components)

# Examples
```jldoctest ; setup = :(using GeometricAlgebra: promote_grades)
julia> promote_grades(4, 0:4, 2, 7)
0:4

julia> promote_grades(4, 0, 2) # even multivectors are worth representing specifically
0:2:4

julia> promote_grades(4, 0, 3) # not worth having a specific type for grades (0, 3) in 4 dims
0:4
```
"""
function promote_grades(dim::Integer, p, q)
	p = promote_grades(dim, p)
	q = promote_grades(dim, q)

	p ⊆ q && return q
	p ⊇ q && return p

	p ⊆ (0, dim) ⊇ q && return (0, dim)
	
	all(iseven, p) && all(iseven, q) && return 0:2:dim
	
	0:dim
end
promote_grades(dim::Integer, p, q, c...) = promote_grades(dim, promote_grades(dim, p, q), c...)

promote_grades(abc::AbstractMultivector{Sig}...) where {Sig} = promote_grades(dimension(Sig), grade.(abc)...)

"""
	resulting_grades(combine, dim, p, q)

Non-zero grade(s) resulting from the application of `combine` on `dim`-dimensional multivectors of grade(s) `p` and `q`.
"""
resulting_grades(combine, dim, P, Q) = promote_grades(dim, (resulting_grades(combine, dim, p::Integer, q::Integer) for p in P, q in Q)...)


# widen grades k to the next narrowest subalgebra:
# scalars ⊂ scalar-pseudoscalar ⊂ even subalgebra ⊂ full algebra
function resulting_grades(::Val{:subalgebra}, dim, k)
	k == 0 && return 0
	k ⊆ (0, dim) && return (0, dim)
	all(iseven, k) ? (0:2:dim) : 0:dim
end



#= Grade projection =#

"""
	componentslice(a, k)

Indices of components of grade(s) `k` in multivector `a`.
"""
function componentslice(a::OrType{<:Multivector}, k::Integer)
	i = findfirst(componentbits(a)) do bits
		count_ones(bits) == k
	end
	isnothing(i) && return 1:0
	range(i, length = ncomponents(dimension(a), k))
end

componentslice(a::OrType{<:Multivector}, k) = reduce(vcat, componentslice.(Ref(a), k); init=Int[])



"""
	a[k]
	getindex(a::Multivector, k)

Get the grade(s) `k` part of a multivector `a` if `k ⊆ grade(a)`.
The components of the resulting `Multivector` are a view into the components of `a`,
so `a` may be modified like `a[k].comps[i] = ...`.
"""
function Base.getindex(a::Multivector, k)
	k = promote_grades(dimension(a), k)
	k ⊆ grade(a) || throw(ArgumentError("""
	attempt to access grade $k part of grade $(grade(a)) Multivector.
	Use `grade(a, k)` to project onto grades which may not exist in `a`.
	"""))
	Multivector{signature(a),k}(view(a.comps, componentslice(a, k)))
end



"""
	grade(a::Multivector{Sig}, k) -> Multivector{Sig,k}

Grade `k` part of a multivector `a`. Returns a zero `k`-vector if `k ∉ grade(a)`.

Multiple grades may be given as a range or tuple. The operators `+` and `-`
may be used as shortcuts for the even and odd parts, respectively.

# Examples
```jldoctest
julia> mv = Multivector{Cl(2,1), 0:4}(1:8)
8-component Multivector{Cl(2,1), 0:4, UnitRange{Int64}}:
 1
 2 v1 + 3 v2 + 4 v3
 5 v12 + 6 v13 + 7 v23
 8 v123

julia> grade(mv, 2)
3-component Multivector{Cl(2,1), 2, UnitRange{Int64}}:
 5 v12
 6 v13
 7 v23

julia> grade(mv, (0, 3))
2-component Multivector{Cl(2,1), (0, 3), SubArray{Int64, 1, UnitRange{Int64}, Tuple{Vector{Int64}}, false}}:
 1
 8 v123

julia> grade(mv, +) # only even grades
4-component Multivector{Cl(2,1), 0:2:2, SubArray{Int64, 1, UnitRange{Int64}, Tuple{Vector{Int64}}, false}}:
 1
 5 v12 + 6 v13 + 7 v23
```

"""
function grade(a::BasisBlade, k)
	k = promote_grades(dimension(a), k)
	T = similar(Multivector{signature(a),k}, a)
	add!(zero(T), a)
end

function grade(a::Multivector{Sig}, k) where {Sig}
	k = promote_grades(dimension(a), k)
	k == grade(a) && return a
	k ⊆ grade(a) && return a[k]

	b = zero(similar(Multivector{Sig,k}, a))
	ak = grade(a) ∩ k
	b.comps[componentslice(b, ak)] = a.comps[componentslice(a, ak)]
	b
end

grade(a::Multivector, ::typeof(+)) = grade(a, 0:2:dimension(a))
grade(a::Multivector, ::typeof(-)) = grade(a, 1:2:dimension(a))

eachgrade(a::Multivector) = (Multivector{signature(a),k}(view(a.comps, componentslice(a, k))) for k in grade(a))
