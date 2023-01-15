```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# GeometricAlgebra

[GeometricAlgebra.jl](https://github.com/jollywatt/GeometricAlgebra.jl) implements flexible types for working with geometric (or Clifford) algebras.

## Quick Start

Construct multivectors by providing a metric signature and grade as type parameters:

```jldoctest
julia> using GeometricAlgebra

julia> u = Multivector{3,1}([1, -1, 0]) # 3D Euclidean vector
3-component Multivector{3, 1, Vector{Int64}}:
  1 v1
 -1 v2
  0 v3
```

Non-euclidean metric signatures may be specified:

```jldoctest
julia> v = Multivector{(-1,1,1,1),2}(1:6) # Lorentzian bivector
6-component Multivector{⟨-+++⟩, 2, UnitRange{Int64}}:
 1 v12
 2 v13
 3 v23
 4 v14
 5 v24
 6 v34

julia> exp(v)
8-component Multivector{⟨-+++⟩, 0:2:4, MVector{8, Float64}}:
 1.18046
 0.818185 v12 + -0.141944 v13 + 0.153208 v23 + 1.076 v14 + 1.16194 v24 + 1.03866 v34
 0.999268 v1234
```

Notice that this bivector exponential has grades `0:2:4`.
The grade parameter `K` of a `Multivector{Sig,K}` can be a single integer
(for homogeneous multivectors) or a collection of grades.
Storing the even part of `exp(v)` (8 components) is more efficient than the full 4D multivector (16 components).

### Working with an orthonormal basis

You may also obtain an orthonormal basis for a metric signature:

```jldoctest
julia> v = basis(3)
3-element Vector{BasisBlade{3, 1, Int64}}:
 v1
 v2
 v3
```

Macros are provided for interactive use:

```jldoctest
julia> @basis "+---"
[ Info: Defined basis blades v1, v2, v3, v4, v12, v13, v23, v14, v24, v34, v123, v124, v134, v234, v1234, I in Main

julia> @basis (t = +1, x = -1) allperms=true
[ Info: Defined basis blades t, x, tx, xt, I in Main
```
