# GeometricAlgebra.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jollywatt.github.io/GeometricAlgebra.jl/dev/)
![Build Status](https://github.com/Jollywatt/GeometricAlgebra.jl/actions/workflows/CI.yml/badge.svg)
[![Coverage](https://codecov.io/gh/jollywatt/GeometricAlgebra.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jollywatt/GeometricAlgebra.jl)

Yet another Julia package for working with geometric (or Clifford) algebras.

## Quick Start

Construct multivectors by providing the metric signature and grade as type parameters:

```julia
julia> u = Multivector{(1,1,1),1}([1, -1, 0])
3-component Multivector{⟨+++⟩, 1, Vector{Int64}}:
  1 v1
 -1 v2
  0 v3

julia> v = Multivector{(1,1,1),2}(1:3)
3-component Multivector{⟨+++⟩, 2, UnitRange{Int64}}:
 1 v12
 2 v13
 3 v23

julia> u*v + π
8-component MixedMultivector{⟨+++⟩, Vector{Float64}}:
 3.14159
 1.0 v1 + 1.0 v2 + -1.0 v3
 5.0 v123
```

You may also obtain an orthonormal basis for a metric signature:

```julia
julia> v = basis((-1,+1,+1,+1))
4-element Vector{Blade{⟨-+++⟩, 1, Int64}}:
 1v1
 1v2
 1v3
 1v4

julia> exp(10000*2π*v[3]v[4])
16-component MixedMultivector{⟨-+++⟩, Vector{Float64}}:
 1.0
 -9.71365e-13 v34
```

Macros are provided for interactive use:

```julia
julia> @basis "+---"
[ Info: Defined basis blades v, v1, v2, v3, v4, v12, v13, v14, v23, v24, v34, v123, v124, v134, v234, v1234

julia> v1234^2
Blade{⟨+---⟩, 0, Int64}:
 -1
```


## Design


There are three concrete types for representing elements in a geometric algebra, arranged in the following type hierarchy:

```
                   AbstractMultivector{Sig}
                     /                  \
   HomogeneousMultivector{Sig,K}    MixedMultivector{Sig,S}
       /                \                             
Blade{Sig,K,T}    Multivector{Sig,K,S}                
                                                   
                  ╰───── CompositeMultivector{Sig,S} ─────╯
```

- `Blade`: a scalar multiple of a wedge product of orthogonal basis vectors.
- `Multivector`: a homogeneous multivector; a sum of same-grade blades.
- `MixedMultivector`: an inhomogeneous multivector. All elements in a geometric
   algebra can be represented as this type (though not most efficiently).

*Note:*
	The mathematical definition of a ``k``-blade is the wedge product
	of ``k`` _vectors_, not necessarily basis vectors. Thus, not all
	``k``-blades are representable as a `Blade`, but are always representable
	as a sum of `Blade`s, or a `Multivector`.

These types have up to three of type parameters:

- `Sig`: The metric signature which defines the geometric algebra. This can be any
   all-bits value which satisfies the metric signature interface.
- `T`: The numerical type of the coefficient of a `Blade`.
- `K`: An `Int` specifying the grade of a `HomogeneousMultivector`.
- `S`: The storage type of the components of a `CompositeMultivector`. This is
   assumed to be mutable, and is usually a subtype of `Vector`, `MVector` or `SparseVector`.

