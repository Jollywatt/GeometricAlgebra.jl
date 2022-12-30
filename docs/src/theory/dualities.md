```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Automorphisms and Dualities


## Fundamental Automorphisms

Generally, operations like complex conjugation ``\overline{AB} = \bar{A}\bar{B}`` or matrix transposition ``(AB)^⊺ = B^⊺A^⊺`` are useful because they preserve or reverse multiplication. (These are called [automorphisms](https://en.wikipedia.org/wiki/Automorphism) and [antiautomorphisms](https://en.wikipedia.org/wiki/Antihomomorphism) respectively.)

Geometric algebras possess some important automorphisms: _reversion_ ``\tilde{A}`` and _grade involution_.

### Reversion

`~A` or [`reversion(a)`](@ref)

Reversion ``\tilde{A}`` is defined on multivectors ``A`` and ``B`` by the property
```math
\tilde{}\,(AB) = \tilde{B}\tilde{A}
```
and by ``\tilde{𝒖} = 𝒖`` for vectors.
Computing the reversion looks like reversing the order of the geometric product:
```@repl ga
~(v1*v2*v3) == v3*v2*v1
```

Swapping orthogonal basis vectors ``\e_i\e_j ↦ \e_j\e_i = -\e_i\e_j`` introduces an overall factor of ``-1``, and it takes ``\binom{k}{2} = \frac{k(k - 1)}{2}`` swaps to reverse ``k`` many basis vectors.
Thus, the reversion of a homogeneous ``k``-vector ``A_k`` is given by
```math
\tilde{A_k} = (-1)^{k(k - 1)/2} A_k
```
but for inhomogeneous multivectors, reversion is not always an overall change in sign.


Grade | Reversion sign
:----:|:-----:
``k`` | ``(-1)^{k(k - 1)/2}``
``0`` | ``+1``
``1`` | ``+1``
``2`` | ``-1``
``3`` | ``-1``
``4`` | ``+1``
``5`` | ``+1``
``6`` | ``-1``
``7`` | ``-1``
``⋮`` | ``⋮``

### [Grade involution](@id involution)

[`involution(A)`](@ref)

Grade involution, sometimes denoted ``A^\star``, is the operation of reflecting space through the origin, so that vectors are sent to their negative, ``𝒖 ↦ -𝒖``. Involution is required to satisfy
```math
\mathsf{involution}(AB) = \mathsf{involution}(A)\mathsf{involution}(B)
```
which means a ``k``-blade of the form ``B = \e_1 ∧ \cdots ∧ \e_k`` gets sent to ``(-\e_1) ∧ \cdots ∧ (-\e_k) = (-1)^k B``.
By linearity, for any ``k``-vector ``A`` we have
```math
\mathsf{involution}(A_k) = (-1)^k A_k
```
but for inhomogeneous multivectors, involution is not always an overall change in sign.

### Clifford conjugation

[`clifford_conj(A)`](@ref)

The composition of reversion and involution ``\tilde{A}^\star`` is also called the Clifford conjugate.



## Pseudoscalars and Dualities

The highest-grade elements, pseudoscalars, play a special role in geometric algebra. The **unit pseudoscalar**
```math
I ≔ \e_1\e_2⋯\e_n
```
is interpreted as an oriented unit volume.

!!! note
	The unit pseudoscalar ``I`` is not be confused with an identity matrix ``𝕀`` or unit imaginary ``i``. Indeed, ``I`` does always commute, and ``I^2 = ±1`` depending on the algebra.

Multiplying by ``I`` sends ``k``-vectors to ``(n - k)``-vectors.