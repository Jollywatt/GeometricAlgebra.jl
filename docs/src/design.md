```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

# Design and Internals

## Multivector Types


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
   all-bits value which satisfies the [metric signature interface](@ref sig-interface).
- `T`: The numerical type of the coefficient of a `Blade`.
- `K`: An `Int` specifying the grade of a `HomogeneousMultivector`.
- `S`: The storage type of the components of a `CompositeMultivector`. This is
   assumed to be mutable, and is usually a subtype of `Vector`, `MVector` or `SparseVector`.


## Metric Signatures

The metric signature type parameter `Sig` defines the dimension of the geometric algebra and the norms of its standard orthonormal basis vectors.
Additionally, it allows various default behaviours to be customised through method definitions which dispatch on `Sig`, as detailed in [the metric signature interface](@ref sig-interface).

By default, the following metric signature types are implemented:
- `Int`, defining a Euclidean metric of that dimension,
- `Tuple`, defining the norms of each basis vector,
- `NamedTuple`, defining basis vector labels as well as norms,
- [`Cl`](@ref), a type resembling the notation ``Cl(p, q, r)`` common in literature.

```jldoctest
julia> @basis 2
[ Info: Defined basis blades v, v1, v2, v12

julia> basis((t=-1, x=1, y=1, z=1)) |> prod
Blade{(t = -1, x = 1, y = 1, z = 1), 4, Int64}:
 1 txyz

```


### [The metric signature interface](@id sig-interface)


The metric signature type parameter may be any `isbits` value satisying the following interface.
As well as defining the geometric algebra, the signature is used to defines basis blade labels, the default array type for multivector components, and so on.

| Required methods | Description |
|:-----------------|:------------|
| `dimension(sig)` | The dimension of the underlying vector space, or number of basis vectors.
| `basis_vector_norm(sig, i)` | The norm of the `i`th basis vector. |

| Optional methods | Description |
|:-----------------|:------------|
| `componentstype(sig, N, T)` | Preferred array type for `CompositeMultivector` components. (Default is `Vector{T}`.)
| `show_signature(io, sig)` | Show the metric signature in a compact form.
| `show_basis_blade(io, sig, indices)` | Print the basis blade with the given indices.

Below is an example of how one might define a geometric algebra with specific behaviours:
```jldoctest
struct DiracGamma end

# define the algebra
GeometricAlgebra.dimension(::DiracGamma) = 4
GeometricAlgebra.basis_vector_norm(::DiracGamma, i) = i > 1 ? -1 : +1

# set the preferred component storage type
using StaticArrays
GeometricAlgebra.componentstype(::DiracGamma, N, T) = MVector{N,T}

# custom labels
GeometricAlgebra.show_basis_blade(io, ::DiracGamma, indices) = print(io, join("γ".*GeometricAlgebra.superscript.(indices)))

basis(DiracGamma())
# output
4-element Vector{Blade{DiracGamma(), 1, Int64}}:
 γ¹
 γ²
 γ³
 γ⁴
```


## Symbolic Algebra and Code Generation

Thanks to the wonderful [`SymbolicUtils`](https://symbolicutils.juliasymbolics.org/) package, the same code originally written for numerical multivectors readily works with symbolic components.
For example, we can compute the product of two vectors symbolically as follows:

```jldoctest
julia> GeometricAlgebra.symbolic_components.([:x, :y], 3)
2-element Vector{Vector{SymbolicUtils.Term{Real, Nothing}}}:
 [x[1], x[2], x[3]]
 [y[1], y[2], y[3]]

julia> Multivector{3,1}.(ans)
2-element Vector{Multivector{3, 1, Vector{SymbolicUtils.Term{Real, Nothing}}}}:
 x[1]v1 + x[2]v2 + x[3]v3
 y[1]v1 + y[2]v2 + y[3]v3

julia> prod(ans)
8-component MixedMultivector{3, Vector{Any}}:
 x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
 x[1]*y[2] - x[2]*y[1] v12 + x[1]*y[3] - x[3]*y[1] v13 + x[2]*y[3] - x[3]*y[2] v23

```

This makes it easy to achieve high performance for basic multivector operations by first performing the calculation symbolically, then converting the resulting expression into unrolled code. (See [`generated_multivector_function`](@ref) for details.)

```julia
julia> u, v = Multivector{4,2}.(eachcol(rand(6,2)));

julia> @time u*v
  0.000009 seconds (2 allocations: 208 bytes)
16-component MixedMultivector{4, Vector{Float64}}:
 -1.31598
 -0.076872 v12 + 0.620294 v13 + 0.344733 v23 + 0.509282 v14 + -0.670917 v24 + -0.619376 v34
 1.34338 v1234

```