# GeometricAlgebra.jl
_A simple, generic multivector implementation for geometric algebra with Julia._

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/Jollywatt/GeometricAlgebra.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jollywatt.github.io/GeometricAlgebra.jl/dev)
[![codecov](https://codecov.io/gh/Jollywatt/GeometricAlgebra.jl/branch/master/graph/badge.svg?token=QP68V6JA6S)](https://codecov.io/gh/Jollywatt/GeometricAlgebra.jl)

The aim of this project is to provide an interface for geometric algebra
which is simple, well documented and easily maintainable.
See also the more mature [Grassmann.jl](https://github.com/chakravala/Grassmann.jl/) and [Multivectors.jl](https://github.com/digitaldomain/Multivectors.jl) packages for more advanced capability.

_This code is very young and in rapid development._

```Julia
julia> x, y, z = basis((1, 1, 1))
3-element Vector{Blade{(1, 1, 1), 1, bits, Int64} where bits}:
 1 v1
 1 v2
 1 v3

julia> x^2 == 1
true

julia> (1 + 2x + 3y)*y
MixedMultivector{(1, 1, 1), Vector{Int64}}:
 3
 0 v1 + 1 v2 + 0 v3
 2 v12 + 0 v13 + 0 v23

julia> ans[2,1]
-2

julia> ~(x*y*z)
Grade-3 Blade{(1, 1, 1), 3, 0x0000000000000007, Int64}:
 -1 v123
```
