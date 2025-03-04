```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

```@setup ga
using BenchmarkTools
using GeometricAlgebra
```

# User Guide

[GeometricAlgebra.jl](https://github.com/jollywatt/GeometricAlgebra.jl) implements a flexible, performant type for multivectors in geometric (or Clifford) algebra.

## Installation

This package is not yet registered, but can be installed directly from source.

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/Jollywatt/GeometricAlgebra.jl", rev="v0.2.1")

julia> using GeometricAlgebra
```

## Quick start

To construct a geometric multivector from a vector of components, provide the _metric signature_ and _grade_ as type parameters.

For example, a 3D Euclidean vector of grade one:
```@repl ga
u = Multivector{3,1}([1, -1, 0])
```
The `Multivector` type is simply a wrapper which associates a geometric algebra and multivector grade(s) to the underlying components vector (accessed with `u.comps`).


Many multivector operations are implemented, including:

- `+`, `-`, `*` ([`geometric_prod`](@ref)), `inv`, `/`
- [`grade`](@ref) for grade projection, [`scalar`](@ref)
- `∧` ([`wedge`](@ref)), `∨` ([`antiwedge`](@ref)), `⋅` ([`inner`](@ref)), `⊙` ([`scalar_prod`](@ref)), `⨼` ([`lcontract`](@ref)), `⨽` ([`rcontract`](@ref))
- `~` ([`reversion`](@ref)), [`involution`](@ref), [`clifford_conj`](@ref)
- [`ldual`](@ref), [`rdual`](@ref), [`hodgedual`](@ref)
- `exp`, `log`, trigonometric functions


Non-euclidean [metric signatures](@ref sig) can be specified, such as `Cl(3,0,1)` for projective geometric algebra (PGA), or a named tuple such as `(t=-1, x=+1)` for custom basis vector names (see [Custom basis display styles](@ref) for more control).

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

## Code generation

The [`@symbolicga`](@ref) macro may be used to generate loop- and allocation-free code from geometric algebra expressions.

```@repl ga
volume(v1, v2, v3) = @symbolicga 3 (v1=1, v2=1, v3=1) scalar(rdual(v1 ∧ v2 ∧ v3))
@btime volume((1, 1, 1), (2, 3, 4), (5, 0, 5))
```

## Custom basis display styles

You can customise both the notation and canonical ordering with which basis blades are displayed to agree with conventions.
This can be done by defining a [`BasisDisplayStyle(dim; kwargs...)`](@ref) and setting a key of [`GeometricAlgebra.BASIS_DISPLAY_STYLES`](@ref).

Some example styles for `v12 + 2v13` are:

| Notation | Display style | `BasisDisplayStyle` keyword arguments
|:--------:|:--------------|:-----------------
| ``𝐞_{12} + 2𝐞_{13}`` | `𝐞12 + 2𝐞13` | `prefix="𝐞"`
| ``γ^0γ^1 + 2γ^0γ^2`` | `γ⁰γ¹ + 2γ⁰γ²` | `prefix="γ", sep="", indices="⁰¹²³"`
| ``\mathrm{d}x ∧ \mathrm{d}y - 2 \mathrm{d}z ∧ \mathrm{d}x`` | `dx ∧ dy - 2 dz ∧ dx` | `prefix="d", sep=" ∧ ", indices="xyz"`

The last style additionally uses a **custom basis blade ordering**.

### Custom basis blade orderings

By default, multivectors are _displayed_ with their components the same sign and in the same order as they are stored.
Internally, the basis vectors in a blade are encoded in binary (see [`bits_to_indices`](@ref) and [`indices_to_bits`](@ref).) For example:
```@repl ga
BasisBlade{4}(42, 0b1101)
```
Multivector components are stored in order of grade, then by binary value. For example, the components of a full 3D multivector correspond to the basis blades:
```@repl ga
GeometricAlgebra.componentbits(Multivector{3,0:3}) .|> UInt8 .|> bitstring
```

However, it can sometimes be convenient to adopt different conventions.
For example, in 3D, it is common to see “cyclical” basis bivectors, ``(𝐯_{23}, 𝐯_{31}, 𝐯_{12})``, as these are the respective duals of ``(𝐯_1, 𝐯_2, 𝐯_3)``.
This style can be achieved as follows:
```@repl ga
style = BasisDisplayStyle(
	3, # dimension of algebra
	Dict(2 => [[2,3], [3,1], [1,2]]), # indices and order of basis bivectors
	prefix="𝐞"
)
GeometricAlgebra.BASIS_DISPLAY_STYLES[3] = style;
```
The second argument of `BasisDisplayStyle` defines the complete list of basis blade indices for each grade.

With this style, it is easier to “read off” the duals of a 3D (bi)vector:
```@repl ga
u = Multivector{3,1}(rand(1:100, 3))
ldual(u)
```
Notice how the vector and dual vector components align nicely.
To recover the default style:
```@repl ga
delete!(GeometricAlgebra.BASIS_DISPLAY_STYLES, 3)
ldual(u)
```

!!! note
	`BasisDisplayStyle` does not affect how components are stored internally. Bear this in mind when accessing the components field of a `Multivector` when using a style with custom ordering.