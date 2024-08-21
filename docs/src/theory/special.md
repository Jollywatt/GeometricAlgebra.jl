```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Inverses, Roots and Logarithms

In general, finding the inverse ``A^{-1}``, square root ``\sqrt{A}`` or logarithm ``\log A`` of a general multivector ``A`` is difficult. However, for certain cases, explicit formulae exist.

## Multivector inverses

Any multivector ``A`` has either no inverse or exactly one inverse ``A^{-1}`` such that ``AA^{-1} = A^{-1}A = 1``.



### Explicit formulae for multivector inverses


For any metric in up to five dimensions, explicit formulae exist for the inverse of a multivector ``A``.
The implementation used in `GeometricAlgebra.jl` is mainly based on [Hitzer2017](@cite) and is described here.

For a multivector ``A ∈ Cl(ℝ^d, ·)`` with metric ``·`` in ``d`` dimensions, let:

- ``Ā`` be the [Clifford conjugate](@ref "Clifford conjugation")
- ``Â`` be the [involute](@ref involution)
- ``Ã`` be the [reverse](@ref "Reversion")
- ``[A]_K`` denote the negation of grades ``k ∈ K``, i.e.,
```math
[A]_K = \sum_{k=0}^d ⟨A⟩_k · \begin{cases}
	-1 & \text{if } k ∈ K \\
	+1 & \text{otherwise}
.\end{cases}
```

| Special case | Formula |
|:-------------|:--------|
| ``A^2 ∈ ℝ`` | ``A^{-1} = \frac{A}{A^2}`` 
| ``d = 3`` | ``A^{-1} = \frac{ĀÂÃ}{AĀÂÃ}``
| ``d = 4`` | ``A^{-1} = \frac{B}{AB}, B = Ā[AĀ]_{3,4}``
| ``d = 5`` | ``A^{-1} = \frac{B}{AB}, B = ĀÂÃ[AĀÂÃ]_{1,4}``


### Faddeev--LeVerrier inverse algorithm

The [Faddeev--LeVerrier inverse algorithm](https://en.wikipedia.org/wiki/Faddeev–LeVerrier_algorithm) may be used to find the inverse of an ``n \times n`` matrix ``A`` in exactly ``n`` steps (with one matrix multiplication per step).
This algorithm can also be used to invert elements in a geometric algebra by considering a linear representation of multivectors as ``n \times n`` matrices. [Dimiter2024](@cite)
This method is implemented in [`inv_flv_method`](@ref).

Crucially, we do not actually need to find a linear representation to use this method: we only need to know that one exists for a given dimension ``n``.
Assume the dimension of the smallest representation is ``\text{repsize}(A)``.

**Algorithm pseudocode**
1. **given** a multivector ``A``
1. ``n := \text{repsize}(A)``
1. ``c_n := 1``
1. ``N \leftarrow 1``
1. **for** ``k \in (n - 1, n - 2, ..., 1, 0)``
1.   ``N \leftarrow N + c_{k + 1} I``
1.   ``c_k := \frac{n}{k - n} A \odot N``
1. **return** ``A^{-1} = -N/c_0``

The inverse exists if and only if ``c_0 \ne 0``.

If ``A`` is a ``d``-dimensional multivector, then:
```math
\text{repsize}(A) = \begin{cases}
	1 & A = ⟨A⟩_0 \\
	2 & \text{$A$ is a pseudoscalar or (pseudo)vector} \\
	2^{\lfloor d/2 \rfloor} & \text{$A$ is even} \\
	2^{\lceil d/2 \rceil} & \text{in general} \\
\end{cases}
```
In some cases, the algorithm works in fewer steps than given above, particularly for multivectors of particular grades.
However, I have not figured out how to determine the true minimum number of steps.
(For example, by inspection, inverting a ``5``-dimensional ``3``-vector requires only ``n = 4`` steps, not ``8`` as suggested by this prescription.)
It's an interesting problem!


## Formulae for multivector square roots

| Special case | Formula |
|:-------------|:--------|
| ``A^2 ∈ ℝ, A^2 < 0`` | ``\sqrt{A} = \frac{A + λ}{\sqrt{2λ}}, λ = \sqrt{-A^2}`` 
| ``A^2 ∈ ℝ, A^2 > 0, I^2 = -1, AI = IA`` | ``\sqrt{A} = \frac{A + Iλ}{(1 + I)\sqrt{λ}}, λ = \sqrt{A^2}`` 

