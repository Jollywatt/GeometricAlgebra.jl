```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Rotors

Rotors are multivectors which describe proper rotations.
The rotors formalism provides an extremely uniform, elegant and efficient description of rotations.
For motivation, it is helpful to consider rotations in simpler algebras.

In the complex plane, a complex number ``z ∈ ℂ`` is rotated about the origin with the mapping ``z ↦ e^{iθ}z``.
Similarly, a [quaternion](https://juliageometry.github.io/Quaternions.jl/) ``q ∈ ℍ`` is rotated about an axis ``n`` with ``q ↦ e^{n/2}qe^{-n/2}``.
Indeed, the double-sided transformation law
```math
A ↦ e^{B/2}Ae^{-B/2}
```
is general to both ``ℂ`` and ``ℍ``, and in fact applies to geometric algebras of any dimension.
Specifically, the [**rotor**](https://en.wikipedia.org/wiki/Rotor_(mathematics)) ``R = e^{B/2}`` describes a rotation in the plane spanned by the bivector ``B`` by an angle described by its magnitude.

## Reflections and orthogonal transformations

To see how the double-sided transformation law arises, note that any orthogonal transformation in $n$ dimensions is the same as at most $n$ reflections (the [Cartan–Dieudonné theorem](https://en.wikipedia.org/wiki/Cartan%E2%80%93Dieudonn%C3%A9_theorem)).
A reflection is described in a geometric algebra by conjugation with an invertible vector.
For instance, the linear map
```math
	A ↦ -𝒗A𝒗^{-1}
```
reflects the multivector $A$ along the vector $𝒗$ (or across the hyperplane with normal $𝒗$).[^1]

[^1]: To see this, consider the case where $A$ is a vector parallel or orthogonal to $𝒗$.

!!! note
	Scaling each $𝒗$ by a non-zero scalar $λ$ does not affect the reflection.
	Therefore, a direct correspondence exists between reflections and _normalised_ vectors $𝒗̂^2 = ±1$, modulo overall sign (since $𝒗̂$ and $-𝒗̂$ describe the same reflection).

By composing reflections as above, we can obtain any orthogonal transformation, acting as
```math
	A ↦ ±RAR^{-1}
```
for some $R = 𝒗_1𝒗_2⋯𝒗_k$.
The overall sign is positive for an even number of reflections (giving a proper rotation), and negative for an odd number.

Without loss of generality, we may use _normalised_ vectors, so that the inverse is
```math
	R^{-1} = 𝒗̂_k^{-1}\cdots 𝒗̂_2^{-1}𝒗̂_1^{-1} = ±\tilde{R}
```
since $𝒗̂^{-1} = ±𝒗̂$.
Hence, an orthogonal transformation is described by
```math
	A ↦ ±RA\tilde{R}
```
where ``R`` is satisfies ``R^{-1} = ±\tilde{R}``.


## Rotor groups

All such multivectors satisfying $R^{-1} = ±\tilde{R}$ taken together form a [group](https://en.wikipedia.org/wiki/Group_(mathematics)) under the geometric product.
This is called the [pin group](https://en.wikipedia.org/wiki/Pin_group):
```math
	\mathsf{Pin}(p, q) ≔ \big\{ R ∈ Cl(p, q) \mid R\tilde{R} = ±1 \big\}
```
There are two “pinors” for every orthogonal transformation, namely $+R$ and $-R$.
Thus, the pin group forms a [double cover](https://en.wikipedia.org/wiki/Covering_space) of the orthogonal group $\mathsf{O}(p,q)$.

Furthermore, the even-grade elements of $\mathsf{Pin}(p,q)$ form a subgroup, called the _spin_ group:
```math
	\mathsf{Spin}(p, q) ≔ \big\{ R ∈ Cl_+(p, q) \mid R\tilde{R} = ±1 \big\}
```
The spin group, in turn, forms a double cover of the special orthogonal group $\mathsf{SO}(p, q)$.

Finally, the additional requirement that $R\tilde{R} = 1$ defines the restricted spinor group, or the **rotor group**:
```math
	\mathsf{Spin}^+(p, q) ≔ \big\{ R ∈ Cl_+(p, q) \mid R\tilde{R} = 1 \big\}
```
The rotor group is a double cover of the restricted special orthogonal group $\mathsf{SO}^+(p, q)$, which is the identity-connected part of $\mathsf{SO}(p, q)$.

The takeaway is that any orthogonal transformation, including reflections, rotations, and combinations of both, can be described within geometric algebra with rotors, no matter the kind of multivector being transformed, and independent of the dimension or signature of the algebra.
In particular, proper rotations are described by **rotors**, or even multivectors satisfying ``R\tilde{R} = 1``.

## The bivector subalgebra

Bivectors play a special role as the _generators_ of rotors.
Because the even subalgebra is closed under the geometric product, the exponential
```math
	e^B = 1 + B + B^2/2 + ⋯ ∈ \mathsf{Spin}^+
```
of a bivector $B$ is always an even multivector, and the reverse ``\tilde{}\,(e^B) = e^{-B}`` is the inverse.
Therefore, ``e^B ∈ \mathsf{Spin}^+`` is a rotor; and indeed, any rotor ``R ∈ \mathsf{Spin}^+`` is of the form
```math
R = e^B
```
for some bivector ``B``.

!!! note "Bivector Lie algebra"
	Formally, bivectors form a Lie algebra under the commutator product $A × B ≔ \frac12(AB - BA)$.
	Indeed, this demonstrates a Lie group–Lie algebra correspondence between the rotor group $\mathsf{Spin}^+$ and bivectors equipped with $×$.

