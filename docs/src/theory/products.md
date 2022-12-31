```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Wedge, Inner and Contraction Products

The geometric product is the fundamental operation in geometric algebra. Together with grade projection, various other “graded products” may be defined by taking different grade projections of the basic product.

To motivate this, consider a ``p``-vector ``A`` and ``q``-vector ``B``. The product contains parts of every grade between the difference ``|p - q|`` and sum ``p + q`` in steps of two:
```math
AB = ⟨AB⟩_{|p - q|} + ⟨AB⟩_{|p - q| + 2} + ⋯ + ⟨AB⟩_{p + q - 2} + ⟨AB⟩_{p + q}
```
Some of these parts are often useful on their own, and so warrant their own name.

These are summarised below, where the grade of the result is shown for each product between a ``p``-vector and ``q``-vector.

Name | Symbol | Resulting grade
----:|:------:|----
scalar | ``\odot`` | ``0``
wedge | ``∧`` | ``p + q``
inner | ``⋅`` | ``\|p - q\|``
left contraction | ``\lcontr`` | ``q - p``
right contraction | ``\rcontr`` | ``p - q``


## Scalar product

[`A ⊙ B`](@ref) or [`scalar_prod(A, B)`](@ref)


The scalar product is simply the scalar part of the geometric product:
```math
A ⊙ B ≔ ⟨AB⟩_0
```
It is also commonly denoted ``A ∗	B`` or ``⟨AB⟩``.


## Wedge product

[`A ∧ B`](@ref) or [`wedge(A, B)`](@ref)

The wedge product may be defined as the highest-grade part of the geometric product.
For example, if ``A`` and ``B`` are multivectors of grade ``p`` and ``q``, respectively, then ``A ∧ B = ⟨AB⟩_{p + q}``.
This definition can be extended linearly to general multivectors as
```math
A ∧ B ≔ \sum_{p,q} \big⟨⟨A⟩_p ⟨B⟩_q\big⟩_{p + q}
```

Between a vector ``𝒖`` and a multivector ``A``, the wedge product may be written as
```math
𝒖 ∧ A = \frac12(𝒖A + A^\star 𝒖)
```
where ``A^\star`` denotes [involution](@ref involution).
Between two vectors, the wedge product is the antisymmetric part of the product:
```math
𝒖 ∧ 𝒗 = \frac12(𝒖𝒗 - 𝒗𝒖)
```

More generally, the wedge product may be written as an antisymmetrised geometric product
```math
𝒖_1 ∧ ⋯ ∧ 𝒖_k = \frac{1}{k!}\sum_{σ ∈ S_k} \operatorname{sign}(σ) 𝒖_{σ(1)} ∧ ⋯ ∧ 𝒖_{σ(k)}
```
where the sum is over all permutations ``σ`` of the indices ``\{1, ..., k\}``.
This shows the connection to [antisymmetric tensors](https://en.wikipedia.org/wiki/Antisymmetric_tensor), but in practice it is much cheaper to compute the highest-grade part of the geometric product!


## Generalised inner product

[`A ⋅ B`](@ref) or [`inner(A, B)`](@ref)

To complement the wedge product, we also define the generalised inner product as the _lowest_-grade part of the geometric product.
For general multivectors, define
```math
A ⋅ B ≔ \sum_{p,q} \big⟨⟨A⟩_p ⟨B⟩_q\big⟩_{|p - q|}
```
Strictly, this should not be confused with the vector inner product, although they are equivalent on vectors.


## Left and right contractions

[`A ⨼ B`](@ref), [`A ⨽ B`](@ref) or [`lcontract(A, B)`](@ref), [`rcontract(A, B)`](@ref)


The left and right contractions are similar to the generalised inner product, except that they do not involve an absolute value (which arguably makes them more ‘uniform’).

```math
\begin{align*}
A \lcontr B &≔ \sum_{p,q} \big⟨⟨A⟩_p ⟨B⟩_q\big⟩_{q - p}
&&\text{(left contraction)}
\\
A \rcontr B &≔ \sum_{p,q} \big⟨⟨A⟩_p ⟨B⟩_q\big⟩_{p - q}
&&\text{(right contraction)}
\end{align*}
```
With a vector ``𝒖``, we have the general formulae
```math
\begin{align*}
𝒖 \lcontr A &= \frac12(𝒖A - A^\star𝒖)
,&
A \rcontr 𝒖 &= \frac12(A𝒖 - 𝒖A^\star)
\end{align*}
```

The link with the inner product can be seen if ``A_p`` is a ``p``-vector and ``B_q`` a ``q``-vector, so that
```math
A_p ⋅ B_q = \begin{cases} A_p \lcontr B_q & p < q \\ A_p \rcontr B_q & p > q \end{cases}
```
and indeed ``A ⋅ B = A \lcontr B = A \rcontr B = A ⊙ B`` if ``p = q``.

Despite their similarities, the contractions are arguably better behaved than the inner product, since identities with the inner product tend to involve grade-based exceptions while identities with contractions tend to hold in full generally. (See [Dorst2002](@cite) for discussion).

For instance, the contractions obey an associativity relation
```math
(A \lcontr B) \rcontr C = A \lcontr (B \rcontr C)
```
for all multivectors ``A, B`` and ``C``, and interact nicely with the wedge product with the identities
```math
\begin{align*}
	(A \lcontr B)I &= A ∧ (BI)
,&	I(A \rcontr B) &= (IA) ∧ B
\\	A \lcontr (B \lcontr C) &= (A ∧ B) \lcontr C
,&	(A \rcontr B) \rcontr C &= A \rcontr (B ∧ C)
\\	A \odot (B \lcontr C) &= (A ∧ B) \odot C
,&	(A \rcontr B) \odot C &= A \odot (B ∧ C)
\end{align*}
```
where ``I`` is the unit pseudoscalar.[^1]

[^1]: Section 3.3, [Wilson2022](@cite)

