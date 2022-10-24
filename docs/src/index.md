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
julia> using GeometricAlgebra

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

```jldoctest
julia> v = basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 v1
 v2
 v3

julia> exp(10000*2π*v[2]v[3])
8-component Multivector{3, Vector{Float64}}:
 1.0
 -9.71365e-13 v23
```

Macros are provided for interactive use:

```jldoctest
julia> @basis "+---"
[ Info: Defined basis blades v, v1, v2, v3, v4, v12, v13, v14, v23, v24, v34, v123, v124, v134, v234, v1234

julia> @basisall (t = +1, x = -1)
[ Info: Defined basis blades t, x, tx, xt
```
