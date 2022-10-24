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
        HomogeneousMultivector{Sig,K}    Multivector{Sig,S}
            /                \                             
BasisBlade{Sig,K,T}    KVector{Sig,K,S}                
                                                        
                       ╰───── CompositeMultivector{Sig,S} ─────╯
```

- `BasisBlade`: a scalar multiple of a wedge product of orthogonal basis vectors.
- `KVector`: a ``k``-vector or homogeneous multivector; a sum of same-grade blades.
- `Multivector`: a general multivector. All elements in a geometric
   algebra can be represented as this type (though not most efficiently).

These types have up to three of type parameters:

- `Sig`: The metric signature which defines the geometric algebra. This can be any
   all-bits value which satisfies the [metric signature interface](@ref sig-interface).
- `T`: The numerical type of the coefficient of a `BasisBlade`.
- `K`: An `Int` specifying the grade of a `HomogeneousMultivector`.
- `S`: The storage type of the components of a `CompositeMultivector`, usually an `AbstractVector` subtype.


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
BasisBlade{(t = -1, x = 1, y = 1, z = 1), 4, Int64}:
 1 txyz

```


### [The metric signature interface](@id sig-interface)


The metric signature type parameter may be any `isbits` value satisying the following interface.
As well as defining the geometric algebra, the signature is used to specify basis blade labels, the default array type for multivector components, and other metadata.

| Required methods | Description |
|:-----------------|:------------|
| `dimension(sig)` | The dimension of the underlying vector space, or number of basis vectors.
| `basis_vector_norm(sig, i)` | The norm of the `i`th basis vector. |

| Optional methods | Description |
|:-----------------|:------------|
| `show_signature(io, sig)` | Show the metric signature in a compact human-readable form.
| `show_basis_blade(io, sig, indices)` | Print a basis blade with the given indices (e.g., `v12` or `𝒆₁∧𝒆₂`).
| `componentstype(sig, N, T)` | Preferred array type for `CompositeMultivector` components. (Default is `Vector{T}` in low dimensions and `SparseVector{T}` otherwise.)
| `symbolic_optim(sig)` | Whether to use symbolic code generation to optimise multivector products. (Default is true for low dimensions.)


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
function GeometricAlgebra.show_basis_blade(io, ::DiracGamma, indices)
	print(io, join("γ".*GeometricAlgebra.superscript.(indices .- 1)))
end

basis(DiracGamma())
# output
4-element Vector{BasisBlade{DiracGamma(), 1, Int64}}:
 γ⁰
 γ¹
 γ²
 γ³
```


## Symbolic Algebra and Code Generation

Thanks to the wonderful [`SymbolicUtils`](https://symbolicutils.juliasymbolics.org/) package, the same code originally written for numerical multivectors readily works with symbolic components.
For example, we can compute the product of two vectors symbolically as follows:

```jldoctest
julia> GeometricAlgebra.symbolic_components.([:x, :y], 3)
2-element Vector{Vector{SymbolicUtils.Term{Real, Nothing}}}:
 [x[1], x[2], x[3]]
 [y[1], y[2], y[3]]

julia> KVector{3,1}.(ans)
2-element Vector{KVector{3, 1, Vector{SymbolicUtils.Term{Real, Nothing}}}}:
 x[1]v1 + x[2]v2 + x[3]v3
 y[1]v1 + y[2]v2 + y[3]v3

julia> prod(ans)
8-component Multivector{3, Vector{Any}}:
 x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
 x[1]*y[2] - x[2]*y[1] v12 + x[1]*y[3] - x[3]*y[1] v13 + x[2]*y[3] - x[3]*y[2] v23

```

This makes it easy to optimize multivector operations by first performing the general calculation symbolically, then converting the resulting expression into unrolled code.
 (See [`symbolic_optim()`](@ref) for details.)

By default, symbolic code generation is used for most products in up to eight dimensions (above which general algebraic expressions become unwieldy). This can be changed a per-algebra basis by defining methods for [`use_symbolic_optim()`](@ref).
