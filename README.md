# GeometricAlgebra.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jollywatt.github.io/GeometricAlgebra.jl/dev/)
![Build Status](https://github.com/Jollywatt/GeometricAlgebra.jl/actions/workflows/CI.yml/badge.svg)
[![Coverage](https://codecov.io/gh/jollywatt/GeometricAlgebra.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jollywatt/GeometricAlgebra.jl)

Yet another Julia package for working with geometric (or Clifford) algebras.

## Quick Start

Construct multivectors by providing the metric signature and grade as type parameters:

```julia
julia> u = KVector([1, -1, 0]) # 3D Euclidean vector
3-component KVector{3, 1, Vector{Int64}}:
  1 v1
 -1 v2
  0 v3

julia> v = KVector{(-1,1,1,1),2}(1:6) # Lorentzian bivector
6-component KVector{⟨-+++⟩, 2, UnitRange{Int64}}:
 1 v12
 2 v13
 3 v23
 4 v14
 5 v24
 6 v34

julia> (v + 1)^2
16-component Multivector{⟨-+++⟩, Vector{Int64}}:
 -48
 2 v12 + 4 v13 + 6 v23 + 8 v14 + 10 v24 + 12 v34
 16 v1234

```

You may also obtain an orthonormal basis for a metric signature:

```julia
julia> v = basis(3)
3-element Vector{Blade{3, 1, Int64}}:
 v1
 v2
 v3

julia> exp(10000*2π*v[2]v[3])
8-component Multivector{3, Vector{Float64}}:
 1.0
 -9.71365e-13 v23
```

Macros are provided for interactive use:

```julia
julia> @basis "+---"
[ Info: Defined basis blades v, v1, v2, v3, v4, v12, v13, v14, v23, v24, v34, v123, v124, v134, v234, v1234

julia> @basisall (t = +1, x = -1)
[ Info: Defined basis blades t, x, tx, xt
```


## Design


There are three concrete types for representing elements in a geometric algebra, arranged in the following type hierarchy:

```
                   AbstractMultivector{Sig}
                     /                  \
   HomogeneousMultivector{Sig,K}    Multivector{Sig,S}
       /               \                             
Blade{Sig,K,T}   KVector{Sig,K,S}                
                                                   
                 ╰─── CompositeMultivector{Sig,S} ───╯
```

- `Blade`: a scalar multiple of a wedge product of orthogonal basis vectors.
- `KVector`: a ``k``-vector or homogeneous multivector; a sum of same-grade blades.
- `Multivector`: a general multivector. All elements in a geometric
   algebra can be represented as this type (though not most efficiently).


## Symbolic Algebra and Code Generation

Thanks to the wonderful [`SymbolicUtils`](https://symbolicutils.juliasymbolics.org/) package, the same code originally written for numerical multivectors readily works with symbolic components.
For example, we can compute the product of two vectors symbolically as follows:

```julia
julia> GeometricAlgebra.symbolic_components.([:x, :y], 3)
2-element Vector{Vector{SymbolicUtils.Term{Real, Nothing}}}:
 [x[1], x[2], x[3]]
 [y[1], y[2], y[3]]

julia> KVector{3,1}.(ans)
2-element Vector{KVector{3, 1, Vector{SymbolicUtils.Term{Real, Nothing}}}}:
 x[1]v1 + x[2]v2 + x[3]v3
 y[1]v1 + y[2]v2 + y[3]v3

julia> prod(ans)
8-component Multivector{3, Vector{Any}}:
 x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
 x[1]*y[2] - x[2]*y[1] v12 + x[1]*y[3] - x[3]*y[1] v13 + x[2]*y[3] - x[3]*y[2] v23

```

This makes it easy to optimize multivector operations by first performing the calculation symbolically, then converting the resulting expression into unrolled code.
By default, symbolic code generation is used for most products in up to eight dimensions (above which general algebraic expressions become unwieldy).

## Similar Packages

This package derives inspiration from many others:

- [ATell-SoundTheory/CliffordAlgebras.jl](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl)
- [chakravala/Grassmann.jl](https://github.com/chakravala/Grassmann.jl)
- [digitaldomain/Multivectors.jl](https://github.com/digitaldomain/Multivectors.jl)
- [MasonProtter/GeometricMatrixAlgebras.jl](https://github.com/MasonProtter/GeometricMatrixAlgebras.jl)
- [serenity4/GeometricAlgebra.jl](https://github.com/serenity4/GeometricAlgebra.jl)
- [velexi-research/GeometricAlgebra.jl](https://github.com/velexi-research/GeometricAlgebra.jl)
- in the future, [JuliaGeometricAlgebra/GeometricAlgebra.jl](https://github.com/JuliaGeometricAlgebra/GeometricAlgebra.jl)
