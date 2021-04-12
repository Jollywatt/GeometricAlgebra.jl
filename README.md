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
julia> t, x, y, z = basis((t=-1, x=1, y=1, z=1)) # create basis vectors for "-+++" spacetime algebra
4-element Vector{Blade{⟨t-,x+,y+,z+⟩, 1, bits, Int64} where bits}:
 1t
 1x
 1y
 1z

julia> (1 + 2x + 3x*y + 4x*y*z)*x
MixedMultivector{⟨t-,x+,y+,z+⟩, Vector{Int64}}:
 2
 1 x + -3 y
 4 yz

julia> grade(ans, 2)
Grade-2 Multivector{⟨t-,x+,y+,z+⟩, 2, Vector{Int64}}:
 0 tx
 0 ty
 0 xy
 0 tz
 0 xz
 4 yz

julia> ans[4,3] # access zy component
-4

julia> R = exp(π/4*x*y) # rotor describing quarter turn in xy plane
MixedMultivector{⟨t-,x+,y+,z+⟩, Vector{Float64}}:
 0.7071067811865476
 0.7071067811865475 xy

julia> ~R*(x + 5z)*R # rotor application R̃aR with reversion ~
MixedMultivector{⟨t-,x+,y+,z+⟩, Vector{Float64}}:
 2.220446049250313e-16 x + 1.0 y + 5.0 z

julia> log(R) ≈ π/4*x*y
true
```
