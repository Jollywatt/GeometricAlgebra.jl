# GeometricAlgebra.jl
_A simple, generic multivector implementation for geometric algebra with Julia._

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/Jollywatt/GeometricAlgebra.jl)
[![Docs](https://img.shields.io/badge/documentation-blue.svg)](https://jollywatt.github.io/GeometricAlgebra.jl/stable)

The aim of this project is to provide an interface for geometric algebra
which is simple, type-agnostic and explicit.
See also the more mature [Grassmann.jl](https://github.com/chakravala/Grassmann.jl/) and [Multivectors.jl](https://github.com/digitaldomain/Multivectors.jl) packages for more advanced capability.

```Julia
julia> x, y, z = basis((1,1,1))
3-element Vector{Blade{⟨+++⟩, Float64, UInt64, 1}}:
 1.0 v₁
 1.0 v₂
 1.0 v₃

julia> 20x - y/2
1-Multivector{⟨+++⟩, Vector{Float64}, 1}
 20.0 v₁
 -0.5 v₂
  0.0 v₃

julia> exp(π/2*x*y)
MixedMultivector{⟨+++⟩, Vector{Float64}}
 6.123233995736766e-17
 1.0 v₁v₂

julia> Multivector{(t=-1, x=1, y=1, z=1)}(2, rand(Int8, 6)) # create random bivector with spacetime metric
2-Multivector{⟨t-,x+,y+,z+⟩, Vector{Int8}, 2}
  36 tx
  51 ty
 -64 xy
   7 tz
   8 xz
 -74 yz

julia> ans[:z,:x] # access the z∧x component by basis label
-8.0

```
Currently supported operations include the geometric `*`, wedge `∧`, dot `⋅` and scalar `∗` products, etc.
Exponentiation `^` and `exp` are supported, but that's about it so far...


Multivector components can be internally stored in (possibly sparse) vectors or in dictionaries:
```Julia
julia> (x + 2y - z).comps
3-element Vector{Float64}:
  1.0
  2.0
 -1.0

julia> t = basis(Multivector{(-1,1,1,1), Dict{UInt8, Int}, 1}, 3)
1-Multivector{⟨-+++⟩, Dict{UInt8, Int64}, 1}
 1 v₃

julia> t.comps
Dict{UInt8, Int64} with 1 entry:
  0x01 => 1

```

Type information is always displayed in full, with the exception of metric signatures, which are pretty-printed unless they appear on their own:

```Julia
julia> x = Multivector{EuclideanSignature(3)}(2, rand(3))
2-Multivector{⟨3+⟩, Vector{Float64}, 2}
 0.5204970431107714 v₁v₂
 0.3520307128412734 v₁v₃
 0.9980489482280224 v₂v₃

julia> signature(x)
⟨3+⟩ = EuclideanSignature(3)
```

Components may be stored in dictionaries with unit blade keys of type `Vector{<:Integer}` or `Vector{<:Symbols}`.
In principle, this flexibility allows for arbitrarily many dimensions.
(Maybe this could find some use? I don't know.)
```Julia

julia> prod(basis.(Blade{EuclideanSignature(100000),Int,Vector{Int},1}, [1, 2, 99999, 100000]))
4-Blade{⟨100000+⟩, Int64, Vector{Int64}, 4}
 1 v₁v₂v₉₉₉₉₉v₁₀₀₀₀₀

julia> ans + 42
MixedMultivector{⟨100000+⟩, Dict{Vector{Int64}, Int64}}
 42
 1 v₁v₂v₉₉₉₉₉v₁₀₀₀₀₀

julia> ans.comps
Dict{Vector{Int64}, Int64} with 2 entries:
  Int64[]               => 42
  [1, 2, 99999, 100000] => 1
```

Operations which do not depend on the dimension of the space still work even if
the metric signature is of undefined dimension.
This enables playground-style interaction, where the dimension need not be specified in advance.

```Julia
julia> a = basis(EuclideanSignature, :a)
1-Blade{EuclideanSignature, Float64, Vector{Symbol}, 1}
 1.0 a

julia> b = basis(EuclideanSignature, :b)
1-Blade{EuclideanSignature, Float64, Vector{Symbol}, 1}
 1.0 b

julia> b*a
2-Blade{EuclideanSignature, Float64, Vector{Symbol}, 2}
 -1.0 ab

julia> (a - b).comps
Dict{Vector{Symbol}, Float64} with 2 entries:
  [:b] => -1.0
  [:a] => 1.0

```
