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

The [Hodge star operator](https://en.wikipedia.org/wiki/Hodge_star_operator) is a metrical duality operation from [exterior algebra](https://en.wikipedia.org/wiki/Exterior_algebra).
For two ``k``-vectors ``A`` and ``B``, the Hodge dual is defined by
```math
A ∧ \mathsf{hodgedual}(B) = ⟨A, B⟩ I
```
where ``⟨A, B⟩`` is the induced inner product on ``k``-vectors.
In the language of geometric algebra, this is
```math
⟨A, B⟩ = A \odot \tilde{B}
```
and the Hodge dual is the same as reversion followed by right-multiplication by the pseudoscalar:
```math
\mathsf{hodgedual}(A) = \tilde{A}I
```

### Comparison with pseudoscalar-duality

For homogeneous multivectors, Hodge duality and pseudoscalar-duality differ only in overall sign.

The square of the Hodge dual is ``\mathsf{hodgedual}^2(A) = (-1)^s(-1)^{k(n - k)} A`` and hence the inverse is
```math
\mathsf{hodgedual}^{-1}(A) = (-1)^s(-1)^{k(n - k)} \mathsf{hodgedual}(A) 
```
where ``s`` is the trace of the metric.[^1] Note that this depends on the grade ``k`` of ``A``.

[^1]: Lemma 6, [Wilson2022](@cite)

By contrast, ``I^2 = ±1`` and hence ``I^{-1} = ±I`` does not depend on the multivector it acts on. (This generally makes pseudoscalar-duality easier to work with algebraically!)