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

!!! note
	The mathematical definition of a ``k``-blade is the wedge product
	of ``k`` _vectors_, not necessarily basis vectors. Thus, not all
	``k``-blades are representable as a `Blade`, but are always representable
	as a sum of `Blade`s, or a `Multivector`.

These types have up to three of type parameters:

- `Sig`: The metric signature which defines the geometric algebra. This can be any
   all-bits value which satisfies the [metric signature interface](@ref sig_interface).
- `T`: The numerical type of the coefficient of a `Blade`.
- `K`: An `Int` specifying the grade of a `HomogeneousMultivector`.
- `S`: The storage type of the components of a `CompositeMultivector`. This is
   assumed to be mutable, and is usually a subtype of `Vector`, `MVector` or `SparseVector`.


## Metric Signatures

The metric signature `Sig` defines the dimension of the geometric algebra and the norms of its standard orthonormal basis vectors.

```jldoctest
julia> sig = (-1,+1,+1,+1)
(-1, 1, 1, 1)

julia> Multivector{sig,1}(1:4)
4-component Multivector{⟨-+++⟩, 1, UnitRange{Int64}}:
 1 v1
 2 v2
 3 v3
 4 v4

julia> signature(ans)
(-1, 1, 1, 1)

```

Additionally, the `Sig` type parameter carries metadata via multiple dispatch: basis blade labels and the default storage type for multivector components may be defined for each metric signature.

```jldoctest
julia> struct DiracGamma end

julia> GeometricAlgebra.dimension(::DiracGamma) = 4

julia> Base.getindex(::DiracGamma, i) = i > 1 ? -1 : +1

julia> GeometricAlgebra.show_basis_blade(io, ::DiracGamma, indices) = print(io, join("γ".*GeometricAlgebra.superscript.(indices)))

julia> prod(basis(DiracGamma()))
Blade{DiracGamma(), 4, Int64}:
 1 γ¹γ²γ³γ⁴

```

### [Metric signature interface](@id sig_interface)

| Required methods | Description |
|:-----------------|:------------|
| `dimension(sig)` | The dimension of the underlying vector space, or number of basis vectors.
| `Base.getindex(sig, i)` | The norm of the `i`th basis vector. |

| Optional methods | Description |
|:-----------------|:------------|
| `Multivector.show_signature(io, sig)` | Show the metric signature in a compact form.
| `Multivector.show_basis_blade(io, sig, indices)` | Print the basis blade with the given indices.