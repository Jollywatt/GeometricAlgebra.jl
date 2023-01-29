```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Dualities


## Pseudoscalar duality

The highest-grade elements, pseudoscalars, play a special role in geometric algebra. The **unit pseudoscalar**
```math
I ≔ 𝐯_1𝐯_2⋯𝐯_n
```
is interpreted as an oriented unit volume.

!!! note
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

For a unit basis blade ``a``, the complements satisfy
```math
\textsf{ldual}(a)\,a = I = a\,\textsf{rdual}(a)
```
They are inverses of each other,
```math
\textsf{ldual}(\textsf{rdual}(a)) = a = \textsf{rdual}(\textsf{ldual}(a))
```
and in odd-dimensional algebras, are identical.

The metric-independence of the left and right duals mean they are appropriate even for degenerate algebras, since a non-zero multivector always has a non-zero dual, even if its Hodge dual is zero.