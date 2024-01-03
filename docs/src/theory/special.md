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



## Formulae for multivector square roots

| Special case | Formula |
|:-------------|:--------|
| ``A^2 ∈ ℝ, A^2 < 0`` | ``\sqrt{A} = \frac{A + λ}{\sqrt{2λ}}, λ = \sqrt{-A^2}`` 
| ``A^2 ∈ ℝ, A^2 > 0, I^2 = -1, AI = IA`` | ``\sqrt{A} = \frac{A + Iλ}{(1 + I)\sqrt{λ}}, λ = \sqrt{A^2}`` 

