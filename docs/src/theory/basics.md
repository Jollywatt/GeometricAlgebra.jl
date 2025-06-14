```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

```@setup ga
using GeometricAlgebra
```

# Introduction

This is a brief overview of the mathematics of geometric algebra, using `GeometricAlgebra.jl` to explain features along the way.

To follow along, install `GeometricAlgebra.jl` and load it with

```julia
julia> using GeometricAlgebra
```


## Definition

Briefly, [geometric algebra](https://en.wikipedia.org/wiki/Geometric_algebra) is what you get when vectors in a space are allowed to be multiplied _freely_,[^1] with the rule that a vector ``𝒖`` multiplied by itself gives its scalar dot product, ``𝒖^2 = 𝒖⋅𝒖``.

What do you get when you multiply two different vectors?
We can simplify expressions using the relation ``𝒖^2 = 𝒖⋅𝒖``, along with the usual associativity ``𝒖(𝒗𝒘) = (𝒖𝒗)𝒘`` and distributivity ``(𝒖 + 𝒗)𝒘 = 𝒖𝒘 + 𝒗𝒘``, but we may still be left with irreducible products, ``𝒖𝒗``.
These are a _new kind of object_, distinct from scalars or vectors. These objects are called **multivectors**.

[^1]: More formally, ``Cl(V, ⋅)`` is the _freest_ unital associative algebra generated by ``V`` satisfying ``𝒖^2 = 𝒖⋅𝒖`` for all ``𝒖 ∈ V``.

What this means becomes clearer after choosing a basis. In geometric algebras, orthonormal vectors ``𝐯_1, ..., 𝐯_n ∈ V`` multiply according to the rules[^2]
- ``𝐯_i^2 = 𝐯_i⋅𝐯_i``,
- ``𝐯_i𝐯_j = -𝐯_j𝐯_i`` where ``i ≠ j``.
The basis vector _norms_ ``𝐯_i⋅𝐯_i`` may be ``+1``, ``-1`` or ``0`` depending on the [inner product](https://en.wikipedia.org/wiki/Inner_product_space) ``⋅`` on ``V``.
The elements ``𝐯_i𝐯_j`` are [**bivectors**](https://en.wikipedia.org/wiki/Bivector) (a kind of oriented plane-element), and higher grade terms ``𝐯_i𝐯_j𝐯_k`` are **trivectors**, etc.

[^2]: The second rule is implied by the first, which is the defining relation ``𝒖^2 = 𝒖⋅𝒖``.

In code, we may obtain an orthonormal basis with [`basis(sig)`](@ref), where `sig` is the [metric signature parameter](@ref sig).
We may also introduce basis variables into the global namespace with the macro [`@basis`](@ref).
```@repl ga
@basis 3
v1^2 # vector square gives its norm, a scalar (or grade 0 blade)
v2*v1 # orthogonal vectors anticommute
```

!!! note "Terminology"
	Geometric algebras are also called (real) [Clifford algebras](https://en.wikipedia.org/wiki/Clifford_algebra) and are denoted ``Cl(V, ⋅)`` or ``Cl(V, Q)``, where ``⋅`` has the [associated quadratic form](https://en.wikipedia.org/wiki/Bilinear_form#Derived_quadratic_form) ``Q``.

	Denote by ``Cl(n)`` the algebra over ``ℝ^n`` where ``⋅`` is the standard Euclidean dot product.
	Write ``Cl(p,q)`` for ``V = ℝ^{p + q}`` where the inner product has ``p`` orthonormal basis vectors with norm ``+1`` and ``q`` with norm ``-1``, and ``Cl(p,q,r)`` if there are an additional ``r`` basis vectors which square to zero.

!!! note "Formal definition"
	A formal construction of the geometric algebra ``Cl(V, ⋅)`` over a vector space ``V`` with inner product ``⋅`` is as the [quotient algebra](https://en.wikipedia.org/wiki/Quotient_ring#For_algebras_over_a_ring) `` Cl(V, ⋅) ≅ V^⊗ \big/ I ``. Here, the tensor algebra ``V^⊗ = ℝ ⊕ V ⊕ (V ⊗ V) ⊕ ⋯`` is reduced by the [ideal](https://en.wikipedia.org/wiki/Ideal_(ring_theory)) ``I`` which identifies all elements of the form ``𝒖 ⊗ 𝒖 - 𝒖⋅𝒖`` where ``𝒖 ∈ V`` with zero.
	The result is a vector space isomorphic to the exterior algebra ``∧V`` but with an algebraic product which mixes degrees.


## Graded structure

Geometric algebras have a _graded structure_: every element has a grade ``k ∈ \{0, 1, 2..., n\}`` or is a sum of elements of differing grades.
As a vector space, ``Cl(V, ⋅)`` is the direct sum of fixed-grade subspaces
```math
Cl(V, ⋅) = ⨁_{k=0}^n ∧^n V = ℝ ⊕ V ⊕ ∧^2 V ⊕ ⋯ ⊕ ∧^n V
```
where ``∧^k V`` is the ``k``th exterior power of the ``n``-dimensional base space ``V``.

Grade zero elements (``∧^0V = ℝ``) are scalars; grade one elements (``∧^1V = V``) are vectors; and grade-``k`` elements (``∧^k V``) are called **``k``-vectors** or **homogeneous multivectors**.
There are no non-zero ``k``-vectors outside the range ``0 ≤ k ≤ n``, so the subspace of ``∧^nV`` contains the highest-grade objects, called **pseudoscalars**.

Grade | Dimension of subspace | Name
:----:|:---------:|:----
``0`` | ``1`` | scalars
``1`` | ``n`` | vectors
``2`` | ``\binom{n}{2}`` | bivectors
``3`` | ``\binom{n}{3}`` | trivectors
``⋮`` | ``⋮`` | ``⋮``
``k`` | ``\binom{n}{k}`` | ``k``-vectors
``⋮`` | ``⋮`` | ``⋮``
``n - 1`` | ``n`` | pseudovectors
``n`` | ``1`` | pseudoscalars
all | ``2^n`` | multivectors

!!! info "Duality"
	Notice that subspaces of grade ``k`` and ``n - k`` have the same dimension.
	Subspaces ``∧^kV`` and their “pseudo”-prefixed counterparts ``∧^{n - k}V`` are associated through [duality](dualities.md).

A ``k``-vector in ``n``-dimensional Euclidean space is represented as a component vector wrapped in the `Multivector{n,k}` type.

```@repl ga
randn(Multivector{3,2})
```

Elements of ``Cl(V, ⋅)`` may consist of parts of differing grade, and when they do they are called **(inhomogeneous) multivectors**.

These are represented as a single contiguous component vector wrapped in a `Multivector` type where the grade parameter `k` is a range of grades.
For example, a general 4D multivector is expressible as `Multivector{4,0:4}`.

If the inner product ``⋅`` is non-Euclidean, the first type parameter `n` is replaced with a [metric signature](@ref sig).

### Grade projection

If ``A ∈ Cl(V, ⋅)`` is a multivector, then denote its grade ``k`` part as
```math
⟨A⟩_k ∈ ∧^kV
\quad\text{(grade projection)}
```
In code, grade projection is achieved with `grade(A, k)`:
```@repl ga
A = Multivector{4,0:4}(rand(-10:10, 2^4))
grade(A, 3)
```

## Blades and multivectors

A general element in a geometric algebra is called a **multivector**, though there is a hierarchy of specific kinds:

```math
\textsf{basis blades} ⊂
\textsf{blades} ⊂
\textsf{$k$-vectors} ⊂
\textsf{multivectors}
```

However, `GeomtricAlgebra.jl` defines only two types: [`BasisBlade`](@ref) and [`Multivector`](@ref).


### Basis blades

Let ``𝐯_1, ..., 𝐯_n`` be the standard orthonormal basis vectors with ``𝐯_i⋅𝐯_i ∈ \{+1,0,-1\}`` and ``𝐯_i⋅𝐯_j = 0`` for ``i ≠ j``.
Products of the form
```math
𝐯_{i_1}𝐯_{i_2}⋯𝐯_{i_k}
```
where ``i_1 < ⋯ < i_k`` are the standard **basis blades**. Reordering the product only introduces overall factors of ``±1``, flipping the **orientation** of the blade.


Scalar multiples of basis blades are represented with the `BasisBlade{Sig,K,T}` type, which encodes indices ``i_1, ..., i_k`` as binary-ones. (See [`bits_to_indices`](@ref) and [`indices_to_bits`](@ref).)
```@repl ga
BasisBlade{4}(42, 0b1011)
```
You can generate all basis blades (of a given grade) for an algebra with [`basis()`](@ref).
```@repl ga
basis(Cl(3,1), 2) # 4-dimensional Lorentzian 2-blades
```

Basis blades are stored with canonically sorted indices. For example, ``𝐯_2𝐯_1`` is represented as ``-𝐯_1𝐯_2``:
```@repl ga
@basis 2 allperms=true
v21
```
While this is always the way blades are represented internally, how they are _displayed_ can be customised through a [`BasisDisplayStyle`](@ref).


#### Multiplication of basis blades

Basis blades are closed under multiplication: a product of basis blades is a basis blade. (This motivates them having their own `BasisBlade` type.)


The geometric product of ``𝐯_I ≡ 𝐯_{i_1}𝐯_{i_2}⋯𝐯_{i_p}`` and ``𝐯_J ≡ 𝐯_{j_1}𝐯_{j_2}⋯𝐯_{j_q}`` can be put into canonical form it two steps:
1. Sort the basis vectors into ``𝐯_{k_1}𝐯_{k_2}⋯𝐯_{k_{p+q}}`` where ``k_1 ≤ ⋯ ≤ k_{p + q}``. Each transposition of adjacent, distinct vectors results in an overall factor of ``-1``. 
2. Simplify adjacent repeated vectors using the fact that ``𝒗^2 = 𝒗⋅𝒗`` is a scalar.

Mathematically, this is
```math
𝐯_I𝐯_J =
	\underbrace{\left(\sum_{k ∈ I ∩ J} 𝐯_k^2\right)}_{\text{factor from squares}}
	\underbrace{\operatorname{sign}(\operatorname{sortperm}(I ⊻ J))}_{\text{sign from swaps}}
	\; 𝐯_{I ⊻ J}
```
where ``I ∩ J`` are the shared indices and ``I ⊻ J`` is the symmetric difference (the indices present in only one blade). In multi-index notation, ``𝐯_{I ⊻ J} ≡ 𝐯_{k_1}𝐯_{k_2}⋯𝐯_{k_m}`` where ``I ⊻ J = \{k_1, ..., k_m\}``. ``\operatorname{sortperm}(I)`` is the permutation that puts the indices ``I`` in ascending order.

You can print a basis blade multiplication table using [`cayleytable`](@ref):
```@repl ga
cayleytable(3)
```


### Blades

While a _basis_ blade is a product of orthogonal _basis_ vectors, the mathematical definition of a **``k``-blade** is a product of ``k`` distinct, mutually orthogonal vectors — or equivalently, a ``k``-blade is a ``∧``-product of ``k`` linearly independent vectors.

Not all blades are representable as basis blades in a given choice of basis. For example, in ``Cl(3)`` the product ``𝐯_1(𝐯_2 + 𝐯_3)`` is a blade (since ``𝐯_1`` and ``𝐯_2 + 𝐯_3`` are orthogonal vectors), but it is a sum of two basis blades, ``𝐯_1𝐯_2 + 𝐯_1𝐯_3``.

Thus, we cannot represent all blades with the `BasisBlade` type.
Instead, we may use the more general `Multivector` type:
```@repl ga
@basis 3
v1*(v2 + v3) # a blade, but not a basis blade
```

### Multivectors

A sum of ``k``-blades is a ``k``-vector, or **homogeneous multivector**.
A sum of blades of differing grade is an **inhomogeneous multivector**.
We may use terminology such as “``\{0,3\}``-vector” to mean an inhomogeneous multivector with both scalar and trivector parts.

The `Multivector{Sig,K}` type represents a `K`-vector, which may be homogeneous (if `K` is an integer) or inhomogeneous (if `K` is a collection, such as `(0, 3)` or `0:3`). 
In all cases, the underlying data is stored in a single vector, while the type parameters define the interpretation.

```@repl ga
@basis 3
4 + 7I # scalar + pseudoscalar; a (0, 3)-vector
ans.comps # underlying components vector
1 + 2v1 + 3v23 # a general multivector
ans.comps
```


!!! note
	(Pseudo)scalars and (pseudo)vectors are always blades, but not all ``k``-vectors are ``k``-blades.
	The simplest example of a homogeneous multivector which isn’t a blade requires four dimensions: the bivector ``𝐯_1𝐯_2 + 𝐯_3𝐯_4`` is not a blade since it cannot be factored as ``𝒖∧𝒗`` for vectors ``𝒖`` and ``𝒗``.


