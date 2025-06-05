```@meta
CurrentModule = GeometricAlgebra.Conformal
DocTestSetup = quote
	using GeometricAlgebra
	using GeometricAlgebra.Conformal
end
```

```@setup cga
using GeometricAlgebra
using GeometricAlgebra.Conformal
```

# Conformal Geometric Algebra


Conformal geometric algebra (CGA) is an elegant projective model of points, lines, circles, planes, spheres, and other objects in $\mathbb{R}^n$ using blades in a geometric algebra over $\mathbb{R}^{n + 2}$.

This package includes a small submodule which defines CGA and a few important operations, including [`standardform`](@ref) for classifying blades and the functions [`ipns`](@ref) and [`opns`](@ref) for obtaining the geometric forms represented by a blade.

You can obtain the standard basis blades for $n$-D CGA using the `CGA{n}` metric signature.
```@repl cga
using GeometricAlgebra.Conformal
basis(CGA{3})
```
Notice that CGA includes two extra basis vectors: `vp` and `vm` squaring to $+1$ and $-1$, respectively.
More generally, any metric signature `Sig` (defining the _base space_) has a conformalisation `CGA{Sig}` with two extra dimensions.

## The representation of points in $\mathbb{R}^n$

In CGA, points $p \in \mathbb{R}^n$ are represented by null vectors (squaring to zero) in $\mathbb{R}^{n + 2}$ given by the map
```math
\begin{align*}
\operatorname{up} &: \mathbb{R}^n \to \mathbb{R}^{n + 2} \\
\operatorname{up}(p) &= n_0 + p + \frac12 p^2 n_\infty
\end{align*}
```
where $n_0$ and $n_\infty$ are null basis vector satisfying $n_0^2 = n_\infty^2 = 0$ and $n_0 \cdot n_\infty = -1$.
In terms of standard basis vectors, they are defined as `n0 = (vm - vp)/2` and `noo = vm + vp` and can be obtained with [`nullbasis`](@ref), [`origin`](@ref) and [`infinity`](@ref).
```@repl cga
n0, noo = Conformal.nullbasis(3);
```

Use the [`up`](@ref) function to lift $1$-vectors from a base space with signature `Sig` into `CGA{Sig}`:
```@repl cga
@basis 3
p = v1 + 2v2 + 3v3
up(p)
ans == n0 + p + 2\p^2*noo # multivectors in Sig and CGA{Sig} are promoted
```
The overall scale of vectors (and blades) in CGA does not matter; it is a _homogeneous_ model.
We can go the other way with the _down_ map, [`dn`](@ref).
```@repl cga
p == dn(up(p)) == dn(100up(p))
```
You may notice that $n_0$ itself is equal to $\operatorname{up}(0)$ and $n_\infty$ is the limit of $\operatorname{up}(x)/\|x\|^2$ as $x$ goes to infinity in any direction.
We interpret $n_0$ as the origin and $n_∞$ as the unique _point at infinity_.

## Standard form of CGA blades

It is useful to represent blades in `CGA{Sig}` in purely terms of objects in the base space `Sig` and the null basis $n_0$ and $n_∞$.

In fact, any blade $X$ in `CGA{Sig}` is of exactly one of the four forms represented by the subtypes of [`CGABlade`](@ref).

| `CGABlade` Subtype | Mathematical form | $X \wedge n_\infty$ | $X \lfloor n_\infty$ |
|---------|:----:|:-------------------:|:--------------------:|
| [`DirectionBlade`](@ref) | $\bm{E} \wedge n_\infty$ | $=0$ | $=0$ |
| [`FlatBlade`](@ref) | $\mathtt{T}_{\bm{p}}[n_0 \wedge \bm{E} \wedge n_\infty]$ | $=0$ | $\ne0$ |
| [`DualFlatBlade`](@ref) | $\mathtt{T}_{\bm{p}}[\bm{E}]$ | $\ne0$ | $=0$ |
| [`RoundBlade`](@ref) | $\mathtt{T}_{\bm{p}}[(n_0 \pm \textstyle{\frac 12} r^2 n_\infty) \wedge \bm{E}]$ | $\ne0$ | $\ne0$ |

Here, $\bm{E}$ is any blade in the base space, $\mathtt{T}_{\bm{p}}$ is the operator which [`translate`](@ref)s by the base space vector $\bm{p}$ and $r$ is a scalar radius parameter.

The [`standardform`](@ref) method converts any blade in `CGA{Sig}` to one of these representations.
```@repl cga
standardform(v12 + 2v1∧noo)
```
This example shows how $𝐯_1𝐯_2 + 2𝐯_2∧n_∞$ can be written as $𝚃_{2𝐯_2}[𝐯_1𝐯_2]$.

## Geometric objects represented by CGA blades

Any CGA blade $X$ can be associated with a subset of the base space $\mathbb{R}^n$ in two ways, providing geometrical interpretations.
We define the _inner_ and _outer product null spaces_[^1]
```math
\begin{align*}
\operatorname{ipns}(X) &\coloneqq \{p \in \mathbb{R}^n \cup \{∞\} \mid \operatorname{up}(p) \mathop\rfloor X = 0 \} \\
\operatorname{opns}(X) &\coloneqq \{p \in \mathbb{R}^n \cup \{∞\} \mid \operatorname{up}(p) \wedge X = 0 \}
\end{align*}
```
which are related by duality (e.g., $\operatorname{ipns} = \operatorname{opns} \circ \operatorname{hodgedual}$) and may be obtained with [`ipns`](@ref) and [`opns`](@ref).

[^1]:
	In an abuse of notation, $\operatorname{up}(∞) = n_∞$.
	More formally, one may define these as subsets of the projective space $Q \subset \mathbf{P}\mathbb{R}^{n + 2}$ of null vectors, $Q = \{[\operatorname{up}(p)] \mid p \in \mathbb{R}^n\} \cup \{[n_\infty]\}$
	where $[\quad]$ is the projective equivalence class.
	Then we have $\operatorname{ipns}(X) = \{x \in Q \mid x \mathop\rfloor X = 0\}$ and similarly for $\operatorname{opns}$.

For a given blade, the inner or outer product null space may be the empty set or an affine $k$-plane (point, line, plane, and so on) or $k$-sphere (point pair, circle, sphere, and so on) in the extended base space $\mathbb{R}^n \cup \{∞\}$.
Geometric objects of these kinds are represented by the [`FlatGeometry`](@ref) and [`RoundGeometry`](@ref) subtypes of [`CGAGeometry`](@ref).


| `CGABlade` Subtype | $\operatorname{ipns}$ | $\operatorname{opns}$ |
|--------------------|:----------------------|:----------------------|
| [`DirectionBlade`](@ref) | [`PointAtInfinity`](@ref) | [`PointAtInfinity`](@ref) |
| [`FlatBlade`](@ref) | [`EmptySet`](@ref) | [`FlatGeometry`](@ref) |
| [`DualFlatBlade`](@ref) | [`FlatGeometry`](@ref) | [`EmptySet`](@ref) |
| [`RoundBlade`](@ref) | [`RoundGeometry`](@ref) | [`RoundGeometry`](@ref) |


For example, below we find that $\operatorname{opns}$ of the outer product of three conformal points is a circle with a centre, direction and real radius.
```@repl cga
p, q, r = up.(randn(Multivector{3,1}, 3));
opns(p∧q∧r)
```
