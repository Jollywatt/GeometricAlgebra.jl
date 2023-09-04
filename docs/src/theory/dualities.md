```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Dualities


## Pseudoscalar duality

The highest-grade elements of geometric algebras are called pseudoscalars. The **unit pseudoscalar**
```math
I ≔ 𝐯_1𝐯_2⋯𝐯_n
```
plays an important role, and may be interpreted as an _oriented unit volume_.

!!! warning
	The unit pseudoscalar ``I`` is not to be confused with the identity matrix or unit imaginary ``i``. Indeed, ``I`` does not always commute, and ``I^2 = ±1`` depending on the algebra.

Multiplying by ``I`` sends ``k``-vectors to ``(n - k)``-vectors.
In odd dimensions, left- and right-multiplication by the unit pseudoscalar is identical:
```math
AI = IA
\quad\text{(in odd dimensions)}
```
for any multivector ``A``.
However, in even dimensions, odd-grade elements anticommute with ``I``.



## Hodge duality

[`hodgedual(A)`](@ref)

The Hodge dual is a combination of reversion and multiplication by ``I``:
```math
\mathsf{hodgedual}(A) = \tilde{A}I
```
This comes from the [Hodge star operator](https://en.wikipedia.org/wiki/Hodge_star_operator) from [exterior algebra](https://en.wikipedia.org/wiki/Exterior_algebra), which is a metrical duality operation implicitly defined on two ``k``-vectors ``A`` and ``B`` by
```math
A ∧ \mathsf{hodgedual}(B) = ⟨A, B⟩ I
```
where ``⟨A, B⟩`` is the induced inner product on ``k``-vectors.
In the language of geometric algebra, this is
```math
⟨A, B⟩ = A \odot \tilde{B} = A \lcontr \tilde{B}
```
and by using the identity ``(A \lcontr \tilde{B})I = A ∧ (\tilde{B}I)`` we have ``A ∧ \mathsf{hodgedual}(B) = A ∧ (\tilde{B}I)`` which shows the equivalence with the explicit definition above.

### Comparison with pseudoscalar-duality

For homogeneous multivectors, Hodge duality and pseudoscalar-duality differ only in overall sign.

The square of the Hodge dual is ``\mathsf{hodgedual}^2(A) = (-1)^s(-1)^{k(n - k)} A`` and hence the inverse is
```math
\mathsf{hodgedual}^{-1}(A) = (-1)^s(-1)^{k(n - k)} \mathsf{hodgedual}(A) 
```
where ``s`` is the trace of the metric.[^1] Note that this depends on the grade ``k`` of ``A``.

[^1]: Lemma 6, [Wilson2022](@cite)

By contrast, ``I^2 = ±1`` and hence ``I^{-1} = ±I`` does not depend on the multivector it acts on. (This generally makes pseudoscalar-duality easier to work with algebraically!)

## Left and right complements

[`ldual`](@ref), [`rdual`](@ref)

The left and right complements are dual operations which do not involve multiplication by the pseudoscalar ``I``, and so are metric independent.
[Some authors](https://rigidgeometricalgebra.org/wiki/index.php?title=Complements) denote the left and right complements by ``\underline{A}`` and ``\overline{A}``, respectively.

For a unit basis blade ``a``, the complements are uniquely defined by
```math
\textsf{ldual}(a)\,a = I = a\,\textsf{rdual}(a)
```
and are extended linearly to act on general multivectors.
They are inverses of each other,
```math
\textsf{ldual}(\textsf{rdual}(a)) = a = \textsf{rdual}(\textsf{ldual}(a))
```
and in odd-dimensional algebras, are identical.

The metric-independence of the left and right duals mean they are useful even in degenerate algebras, since a non-zero multivector always has a non-zero dual even when its Hodge dual is zero.

## Dualized products

It is common to work both with multivectors and their duals in the same context, and hence convenient to define “dualized” versions of certain products.
For example, the regressive product ``∨`` is the dualized [wedge product](@ref wedge) ``∧``, defined by
```math
D(a ∨ b) = D(a) ∧ D(b)
```
where ``D`` is a dual operation such as ``\textsf{ldual}`` or ``\textsf{rdual}`` (each results in an equivalent definition).

In the context of [point-based projective geometric algebra](https://rigidgeometricalgebra.org/), the ``∧``-product of two objects forms the _join_ (e.g., the line joining two points) while the ``∨``-product of two objects forms the _meet_ (e.g., the point where two lines meet) and vice versa in plane-based PGA.

Similar products may be defined by dualizing the inner product and contractions[^2] but these are less common.

[^2]: See the “anti-products” defined in [Lengyel’s wiki](http://projectivegeometricalgebra.org/).