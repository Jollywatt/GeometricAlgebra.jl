<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./docs/src/assets/logo-dark.svg">
  <img alt="logo" width="120" src="./docs/src/assets/logo.svg">
</picture>

# GeometricAlgebra.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jollywatt.github.io/GeometricAlgebra.jl/dev/)
![Build Status](https://github.com/Jollywatt/GeometricAlgebra.jl/actions/workflows/CI.yml/badge.svg)
[![Coverage](https://codecov.io/gh/jollywatt/GeometricAlgebra.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jollywatt/GeometricAlgebra.jl)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/jollywatt/GeometricAlgebra.jl)

A Julia package for working with multivectors from geometric (or Clifford) algebra.

## Installation

This package readily usable, but unregistered. (Mostly because I haven’t thought of a good package name appropriate for the General registry!)

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/Jollywatt/GeometricAlgebra.jl", rev="v0.2.2")
```

## Quick Start

Construct multivectors by providing a metric signature and grade as type parameters:

```julia
julia> using GeometricAlgebra

julia> u = Multivector{3,1}([1, -1, 0]) # 3D Euclidean vector
3-component Multivector{3, 1, Vector{Int64}}:
  1 v1
 -1 v2
  0 v3
```

Non-euclidean metric signatures may be specified:

```julia
julia> v = Multivector{Cl("-+++"),2}(1:6) # Lorentzian bivector
6-component Multivector{Cl("-+++"), 2, UnitRange{Int64}}:
 1 v12
 2 v13
 3 v23
 4 v14
 5 v24
 6 v34

julia> exp(v)
8-component Multivector{Cl("-+++"), 0:2:4, MVector{8, Float64}}:
 1.18046
 0.818185 v12 + -0.141944 v13 + 0.153208 v23 + 1.076 v14 + 1.16194 v24 + 1.03866 v34
 0.999268 v1234
```

Notice that this bivector exponential has grades `0:2:4`.
The grade parameter `K` of a `Multivector{Sig,K}` can be a single integer
(for homogeneous multivectors) or a collection of grades.
A general 4D multivector has grades `0:4`, but an even multivector
may be more efficiently represented with grades `0:2:4`.

You may also obtain an orthonormal basis for a metric signature:

```julia
julia> v = basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 v1
 v2
 v3

julia> exp(10000*2π*v[2]v[3])
4-component Multivector{3, 0:2:2, Vector{Float64}}:
 1.0
 -9.71365e-13 v23
```

Macros are provided for interactive use:

```julia
julia> @basis "+---"
[ Info: Defined basis blades v1, v2, v3, v4, v12, v13, v23, v14, v24, v34, v123, v124, v134, v234, v1234, I in Main

julia> @basis (t = +1, x = -1) allperms=true
[ Info: Defined basis blades t, x, tx, xt, I in Main
```


## Design

There are two concrete types for representing elements in a geometric algebra:

```
         AbstractMultivector{Sig}
            /               \                             
BasisBlade{Sig,K,T}    Multivector{Sig,K,S}
```

- `BasisBlade`: a scalar multiple of a wedge product of orthogonal basis vectors.
- `Multivector`: a homogeneous or inhomogeneous multivector; a sum of basis blades.

Type parameters:

- `Sig`: The metric signature which defines the geometric algebra. This can be any all-bits value which satisfies the metric signature interface.
- `K`: The grade(s) of a multivector. For `BasisBlade`s, this is an integer, but for `Multivector`s, it may be a collection (e.g., `0:3` for a general 3D multivector).
- `T`: The numerical type of the coefficient of a `BasisBlade`.
- `S`: The storage type of the components of a `Multivector`, usually an `AbstractVector` subtype.


## Symbolic Algebra and Code Generation

Thanks to the wonderful [`SymbolicUtils`](https://symbolicutils.juliasymbolics.org/) package, the same code originally written for numerical multivectors readily works with symbolic components.
For example, we can compute the product of two vectors symbolically as follows:

```julia
julia> GeometricAlgebra.MiniCAS.variables.([:x, :y], 3)
2-element Vector{Vector{Any}}:
 [x[1], x[2], x[3]]
 [y[1], y[2], y[3]]

julia> Multivector{3,1}.(ans)
2-element Vector{Multivector{3, 1, Vector{Any}}}:
 x[1]v1 + x[2]v2 + x[3]v3
 y[1]v1 + y[2]v2 + y[3]v3

julia> prod(ans)
4-component Multivector{3, 0:2:2, Vector{Any}}:
 x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
 x[1]*y[2] - x[2]*y[1] v12 + x[1]*y[3] - x[3]*y[1] v13 + x[2]*y[3] - x[3]*y[2] v23

```

This makes it easy to optimize multivector operations by first performing the calculation symbolically, then converting the resulting expression into unrolled code.
By default, symbolic code generation is used for most products in up to eight dimensions (above which general algebraic expressions become unwieldy).

## Similar Packages

This package derives inspiration from many others. Here is a list of Julia implementations (of varying completeness) which I have come across:

- [ATell-SoundTheory/CliffordAlgebras.jl](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl)
- [brainandforce/CliffordNumbers.jl](https://github.com/brainandforce/CliffordNumbers.jl)
- [chakravala/Grassmann.jl](https://github.com/chakravala/Grassmann.jl)
- [digitaldomain/Multivectors.jl](https://github.com/digitaldomain/Multivectors.jl)
- [MasonProtter/GeometricMatrixAlgebras.jl](https://github.com/MasonProtter/GeometricMatrixAlgebras.jl)
- [serenity4/GeometricAlgebra.jl](https://github.com/serenity4/GeometricAlgebra.jl)
- [serenity4/SymbolicGA.jl](https://github.com/serenity4/SymbolicGA.jl)
- [velexi-research/GeometricAlgebra.jl](https://github.com/velexi-research/GeometricAlgebra.jl)
- in the future, [JuliaGeometricAlgebra/GeometricAlgebra.jl](https://github.com/JuliaGeometricAlgebra/GeometricAlgebra.jl)
