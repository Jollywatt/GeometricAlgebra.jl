```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# GeometricAlgebra

[GeometricAlgebra.jl](https://github.com/jollywatt/GeometricAlgebra.jl) implements basic types for working with geometric (or Clifford) algebras.

## Quick Start

Construct multivectors by providing the metric signature and grade as type parameters:

```jldoctest
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

```jldoctest
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

```jldoctest
julia> @basis "+---"
[ Info: Defined basis blades v, v1, v2, v3, v4, v12, v13, v14, v23, v24, v34, v123, v124, v134, v234, v1234

julia> v1234^2
Blade{⟨+---⟩, 0, Int64}:
 -1
```
