```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Background

This is a brief overview of the mathematics of geometric algebra, using `GeometricAlgebra.jl` to explain features along the way.

## Features of geometric algebra

Geometric algebra unifies many concepts in geometry, linear algebra and physics by offering an elegant language for geometric objects and the transformations on them. In particular, it
- provides a uniform formalism for rotations in any number of dimensions;
- efficiently models Euclidean, projective, and conformal geometries;
- unifies the dot product ``⋅`` and wedge product ``∧`` from exterior calculus into the same algebra, creating an _invertible_ geometric product;
- consists of elements with concrete geometric interpretations;
- generalises the complex numbers ``ℂ``, quaternions ``ℍ``, Pauli and Dirac matrix algebras and realises them as algebras over _real_ vector spaces.


## Abstract definition

Geometric algebra is what you get when vectors in a space are allowed to be multiplied freely, subject to:
- associativity, ``𝒖(𝒗𝒘) = (𝒖𝒗)𝒘``; and
- the condition that a vector ``𝒖`` multiplied by itself is the dot product, ``𝒖^2 = 𝒖⋅𝒖``.

Notice that the ingredients for a geometric algebra are a vector space ``V`` and an inner product ``⋅`` on ``V``.

!!! details
	Formally, the geometric algebra ``𝒢(V, ⋅)`` over a vector space ``V`` with the inner product ``⋅`` is the quotient
	```math
	𝒢(V, ⋅) ≅ (ℝ ⊕ V ⊕ (V ⊗ V) ⊕ ⋯) \big/ \{𝒖 ⊗ 𝒖 - 𝒖⋅𝒖 \mid 𝒖 ∈ V\}
	```
	which the free tensor algebra ``ℝ ⊕ V ⊕ (V ⊗ V) ⊕ ⋯`` except with all elements of the form ``𝒖 ⊗ 𝒖 - 𝒖⋅𝒖`` for ``𝒖 ∈ V`` set to zero.

All other features a geometric algebra follow from the definition, even if they are not obvious at first!
If ``𝒗_1, ..., 𝒗_n ∈ V`` are orthonormal basis vectors, then it turns out that in the geometric algebra ``𝒢(V, ⋅)``:
- Elements are generally non-commutative. E.g., ``𝒗_1𝒗_2 = -𝒗_2𝒗_1``.
- There are more elements than just vectors in ``V``: now there are _multivectors_ of grades ``0, ..., n``. E.g., ``𝒗_1𝒗_2`` is a bivector.
- Multiplication can produce multivectors of _mixed grade_. E.g., if ``𝒖 = 2𝒗_1 + 𝒗_2`` is a vector, then ``𝒖𝒗_2 = 1 + 2𝒗_1𝒗_2`` is a ``\{0, 2\}``-multivector with scalar and bivector parts.
- Most elements are invertible. E.g., if ``𝒖`` is the vector above, then ``𝒖^{-1} = \frac25 𝒗_1 + \frac15 𝒗_2`` and ``𝒖𝒖^{-1} = 1``.
- Other products, including the wedge ``∧``, can be defined in terms of the geometric product and _grade projection_.
  E.g., ``𝒖∧𝒗_2 = 2𝒗_1𝒗_2`` takes the highest-grade part of the product ``𝒖𝒗_2 = 1 + 2𝒗_1𝒗_2``.

## Terminology

| Term | Definition | Example |
|------|:-----------|:--------|
| _multivector_ | Any element of a geometric algebra | ``1, 𝒗_1, 1 + 𝒗_1𝒗_2``
| _``k``-vector_ | A homogeneous multivector, containing only grade-``k`` parts | ``𝒘 = 𝒗_1𝒗_2 + 𝒗_3𝒗_4`` is a ``2``-vector or _bivector_
| _vector_ | Usually refers to a ``1``-vector, even though every multivector is a “vector” in the linear algebra sense | ``𝒖 = 2𝒗_1 + 𝒗_2``
| _``k``-blade_ | A ``k``-vector which can be factored as a ``∧``-product of ``k``-many vectors | ``𝒖∧𝒗_3 = 2𝒗_1𝒗_3 + 𝒗_2𝒗_3`` is a ``2``-blade, but ``𝒘`` is not
| _basis blade_ | A ``k``-blade which is a scalar multiple of the product of orthonormal basis vectors | ``1, 2𝒗_1𝒗_3, -𝒗_1𝒗_2𝒗_3``


## Example in ``ℝ^3``

To see what the definition means, we’ll look at the example of 3D space ``V = ℝ^3`` with the usual Euclidean dot product: the 3D geometric algebra ``𝒢(ℝ^3, ⋅)`` or ``𝒢(3)``.

We may get the standard orthonormal basis ``𝒗_1, 𝒗_2, 𝒗_3 ∈ ℝ^n`` using
```jldoctest 3d
julia> v1, v2, v3 = basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 v1
 v2
 v3
```
From the definition ``𝒖^2 = 𝒖⋅𝒖``, we know how to multiply ``𝒗_1`` with itself:
```jldoctest 3d
julia> v1*v1
BasisBlade{3, 0, Int64}:
 1
```
But what about `v1*v2`? The relation ``𝒖^2 = 𝒖⋅𝒖`` doesn’t immediately say how to multiply two different vectors — but everything about geometric product follows from it! Notice that if we apply the relation to ``𝒖 + 𝒗`` instead of ``𝒖``, we find that
```math
\begin{align*}
	(𝒖 + 𝒗)^2 &= (𝒖 + 𝒗)⋅(𝒖 + 𝒗)
\\	𝒖^2 + 𝒖𝒗 + 𝒗𝒖 + 𝒗^2 &= 𝒖⋅𝒖 + 2𝒖⋅𝒗 + 𝒗⋅𝒗
\\	𝒖𝒗 + 𝒗𝒖 &= 2𝒖⋅𝒗
\end{align*}
```
In our case, ``𝒗_1⋅𝒗_2 = 0`` are orthogonal, so this means
```math
𝒗_1𝒗_2 = -𝒗_2𝒗_1
```
which we can check with code:
```jldoctest 3d
julia> v1*v2
BasisBlade{3, 2, Int64}:
 1 v12

julia> v2*v1
BasisBlade{3, 2, Int64}:
 -1 v12
```
These is a new kind of object — they are grade-2 bivectors in 3D space:
```jldoctest 3d
julia> dimension(v1*v2), grade(v1*v2)
(3, 2)
```

In fact, there are ``2^3 = 8`` different linearly independent objects in the geometric algebra ``𝒢(3)``:
```jldoctest 3d
julia> basis(3, grade=:all)
8-element Vector{BasisBlade{3, _A, Int64} where _A}:
 1
 v1
 v2
 v3
 v12
 v13
 v23
 v123
```
We can add them all to the namespace with [`@basis`](@ref):
```jldoctest 3d
julia> @basis 3
[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123
```

### Multivectors