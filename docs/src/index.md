```@setup ga
using GeometricAlgebra
```

```@setup ga3d
using GeometricAlgebra
x, y, z = basis((x=1, y=1, z=1))
```

# GeometricAlgebra.jl


[`GeometricAlgebra.jl`](https://github.com/Jollywatt/GeometricAlgebra.jl) aims to provide a simple interface for [geometric algebra](https://en.wikipedia.org/wiki/Geometric_algebra) (a.k.a. [Clifford algebra](https://en.wikipedia.org/wiki/Clifford_algebra)) in the Julia language.

```@repl
using GeometricAlgebra
x, y, z = basis((1, 1, 1)); # basis for the algebra of physical space
a = x + 2y + 3z # vector
R = exp(π/4*x*y) # rotor
R*a*~R # rotate the vector
```



## Quick start

### Generate basis vectors for algebra of given signature

```@repl ga
basis((1, 1, 1)) # also `basis(3)`

basis((x=-1, y=-1, z=-1)) # labelled anti-Euclidean basis

basis(OffsetSignature((t=-1, x=+1, y=+1, z=+1), 0:3)) |> sum
```

### Perform arithmetic on blades

```@repl ga3d
x*y*z # a grade 3 blade 

x - 3y # a vector, or grade 1 multivector

(1 - 2x*y)/(1 + x*y) # an inhomogeneous multivector of mixed grade

ans^10
```

### Other products and operations

| math | code | description
|------|------|:-----------
|``a∧b``|`a∧b`, `wedge(a, b)`|wedge product
|``a⋅b``|`a⋅b`|interior (or "fat" dot) product
|``a∗b ≡ ⟨ab⟩``|`a∗b == scalar(a*b)`|scalar product (typed `\ast<tab>`)
|``\tilde{a}``|`~a`|reversion
|``a⨼b``|`a⨼b`|left contraction (`\intprod<tab>`)
|``a⨽b``| `a⨽b`|right contraction (`\intprodr<tab>`)
|``⟨a⟩_k``|`grade(a, k)`|grade production

### Component access

```@repl ga
v = basis(3)
a = 1 + 2v[1] + 3v[1]v[2]+ 4v[1]v[2]v[3]
a[2,1] # access the yx = -xy component
a[] == scalar(a) == 1
pseudoscalar(a)
```

For labelled signatures such as `(x=1, y=1, z=1)`, you can access components using the basis labels, too:

```@repl ga
x, y, z = basis((x=1, y=1, z=1))
a = x - y
a[:x]
a[2]
```

Especially in the theory of relativity, convention dictates that components of spacetime 4-vectors are numbered from zero.
This is possible with an [`OffsetSignature`](@ref):

```@repl ga
STA = OffsetSignature((+1,-1,-1,-1), 0:3)
v = basis(STA)
a = rand(UInt8, 4)'v
a[0] # first ("time") component
```
