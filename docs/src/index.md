```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

```@setup ga
using GeometricAlgebra
```

# User Guide

[GeometricAlgebra.jl](https://github.com/jollywatt/GeometricAlgebra.jl) implements flexible types for working with geometric (or Clifford) algebras.

## Installation

This package is not yet registered, but can be installed directly from source.

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/Jollywatt/GeometricAlgebra.jl", rev="v0.2.1")

julia> using GeometricAlgebra
```

## Quick start

To construct a geometric multivector from a vector of components, provide the _metric signature_ and _grade_ as type parameters.

For example, a 3D Euclidean vector of grade one is a `Multivector{3,1}`.

```@repl ga
u = Multivector{3,1}([1, -1, 0])
```
The `Multivector` type is essentially a wrapper which signifies the parent algebra and grade(s) of the underlying components vector (accessed with `u.comps`).


Many multivector operations are implemented, including:

- `+`, `-`, `*` ([`geometric_prod`](@ref)), `inv`, `/`
- [`grade`](@ref) for grade projection, [`scalar`](@ref)
- `∧` ([`wedge`](@ref)), `∨` ([`antiwedge`](@ref)), `⋅` ([`inner`](@ref)), `⊙` ([`scalar_prod`](@ref)), `⨼` ([`lcontract`](@ref)), `⨽` ([`rcontract`](@ref))
- `~` ([`reversion`](@ref)), [`involution`](@ref), [`clifford_conj`](@ref)
- [`ldual`](@ref), [`rdual`](@ref), [`hodgedual`](@ref)
- `exp`, `log`, trigonometric functions


Non-euclidean metric signatures can be specified, such as `Cl(3,0,1)` for projective geometric algebra (PGA), or a named tuple such as `(t=-1, x=+1)` for custom basis vector names (see [`BasisDisplayStyle`](@ref) for more control).

For example, here is a bivector in the spacetime algebra (STA) using the ``({-}{+}{+}{+})`` metric.
```@repl ga
B = Multivector{Cl("-+++"),2}(1:6) # Lorentzian bivector
exp(B)
```
Notice that the grades parameter of this bivector exponential is `0:2:4`, for an even multivector having 8 components (instead of 16 components of a general 4D multivector).


## Working with an orthonormal basis

You can obtain an orthonormal basis for a metric signature with [`basis`](@ref).

```@repl ga
v = basis(3)
```

`BasisBlade`s represent individual components of a `Multivector`.
The macro [`@basis`](@ref) introduces basis blades into the current namespace for interactive use.

```@repl ga
@basis "+---"
1 + 10v23 + 1000/I
@basis (t = +1, x = -1) allperms=true
6tx == -6xt
```

## Custom basis display styles

There are many notations used in geometric algebra literature, so the way basis blades are printed can be customised.
This is done by constructing a [`BasisDisplayStyle`](@ref). E.g., `v12 + 2v13` might be written as

| Notation | Display style | Keyword arguments
|:--------:|:--------------|:-----------------
| ``𝐞_{12} + 2𝐞_{13}`` | `𝐞12 + 2𝐞13` | `prefix="𝐞"`
| ``γ^0γ^1 + 2γ^0γ^2`` | `γ⁰γ¹ + 2γ⁰γ²` | `prefix="γ", sep="", indices=0:3`
| ``\mathrm{d}x ∧ \mathrm{d}y - 2 \mathrm{d}z ∧ \mathrm{d}x`` | `dx ∧ dy - 2 dz ∧ dx` | `prefix="d", sep=" ∧ ", indices="xyz` with a custom ordering

The last style uses a **custom basis blade ordering**.

By default, multivectors are _displayed_ with their components the same sign and in the same order as they are stored.
Internally, basis blades are encoded in binary.
```@repl ga
BasisBlade{4}(42, 0b1101)
```
Multivector components are stored in order of grade, then binary value. For example, the components of a full 3D multivector correspond to basis blades as follows:
```@repl ga
GeometricAlgebra.componentbits(Multivector{3,0:3}) .|> UInt8 .|> bitstring
```

However, it can be convenient to adopt different conventions in certain situations.
For example, in 3D, it is useful to use a “cyclical” basis of bivectors, ``(𝐯_{23}, 𝐯_{31}, 𝐯_{12})``, which are the respective duals of ``(𝐯_1, 𝐯_2, 𝐯_3)``.
This style can be achieved as follows:
```@repl ga
style = BasisDisplayStyle(
  3, # dimension of algebra
  Dict(2 => [[2,3], [3,1], [1,2]]), # indices and order of basis bivectors
  prefix="𝐞"
)
GeometricAlgebra.BASIS_DISPLAY_STYLES[3] = style;
```
Now it is easier to “read off” the duals of a vector or bivector:
```@repl ga
u = Multivector{3,1}(rand(1:100, 3))
ldual(u)
```
Notice how the vector and dual vector components correspond clearly.
To recover the default style:
```@repl ga
delete!(GeometricAlgebra.BASIS_DISPLAY_STYLES, 3)
ldual(u)
```

