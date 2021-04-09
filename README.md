# GeometricAlgebra.jl
_A simple, generic implementation for geometric (Clifford) algebras in Julia._

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/Jollywatt/GeometricAlgebra.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jollywatt.github.io/GeometricAlgebra.jl/dev)
![build](https://api.travis-ci.com/Jollywatt/GeometricAlgebra.jl.svg)
[![codecov](https://codecov.io/gh/Jollywatt/GeometricAlgebra.jl/branch/master/graph/badge.svg?token=QP68V6JA6S)](https://codecov.io/gh/Jollywatt/GeometricAlgebra.jl)

The aim of this project is to provide an interface for geometric algebra
which is simple, well documented and easily maintainable.
See also the more mature [Grassmann.jl](https://github.com/chakravala/Grassmann.jl/) and [Multivectors.jl](https://github.com/digitaldomain/Multivectors.jl) packages for more advanced capability.

_This code is very young and in rapid development._

```Julia
julia> x,y,z = basis((x=1,y=1,z=1))
3-element Vector{Blade{⟨x+,y+,z+⟩, 1, bits, Int64} where bits}:
 1x
 1y
 1z

julia> (1 + 2x + 3x*y + 4x*y*z)x
MixedMultivector{⟨x+,y+,z+⟩, Vector{Int64}}:
 2
 1 x + -3 y + 0 z
 0 xy + 0 xz + 4 yz

julia> grade(ans, 2)
Grade-2 Multivector{⟨x+,y+,z+⟩, 2, Vector{Int64}}:
 0 xy
 0 xz
 4 yz

julia> ans[3,2]
-4

julia> R = exp(π/4*x*y)
MixedMultivector{⟨x+,y+,z+⟩, Vector{Float64}}:
 0.7071067811865476
 0.7071067811865476 xy + 0.0 xz + 0.0 yz

julia> ~R*(x + 5z)*R
MixedMultivector{⟨x+,y+,z+⟩, Vector{Float64}}:
 0.0 x + 1.0000000000000002 y + 5.000000000000001 z

julia> log(R) ≈ π/4*x*y
true
```
