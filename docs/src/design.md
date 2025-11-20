```@meta
CurrentModule = GeometricAlgebra
DocTestSetup = quote
	using GeometricAlgebra
end
```

```@setup ga
using GeometricAlgebra
```

# Design and Internals

## Multivector Types


There are two concrete types for representing elements in a geometric algebra:

```
        AbstractMultivector{Sig,K}
            /               \
BasisBlade{Sig,K,T}    Multivector{Sig,K,S}
```

- [`BasisBlade`](@ref): a scalar multiple of a wedge product of orthogonal basis vectors.
- [`Multivector`](@ref): a homogeneous or inhomogeneous multivector; a sum of basis blades.

Type parameters:

- `Sig`: The metric signature which defines the geometric algebra. This can be any
   all-bits value which satisfies the [metric signature interface](@ref sig-interface).
- `K`: The grade(s) of a multivector. For `BasisBlade`s, this is an integer, but for `Multivector`s, it may be a collection (e.g., `0:3` for a general 3D multivector).
- `T`: The numerical type of the coefficient of a `BasisBlade`.
- `S`: The storage type of the components of a `Multivector`, usually an `AbstractVector` subtype.


## [Metric Signatures](@id sig)

The metric signature type parameter `Sig` defines the dimension of the geometric algebra and the squares of its standard orthonormal basis vectors.
Additionally, it allows various default behaviours to be customised through method definitions which dispatch on `Sig`, as detailed in [the metric signature interface](@ref sig-interface).

By default, the following metric signature types are implemented:
- `Int`, defining a Euclidean metric of that dimension,
- `Tuple`, defining the norms of each basis vector,
- `NamedTuple`, defining basis vector labels as well as norms,
- [`Cl`](@ref), a type resembling the notation ``Cl(p, q, r)`` common in literature.

```jldoctest
julia> @basis 2
[ Info: Defined basis blades v1, v2, v12, I in Main

julia> basis((t=-1, x=1, y=1, z=1)) |> prod
BasisBlade{(t = -1, x = 1, y = 1, z = 1), 4, Int64}:
 1 txyz

julia> sum(basis(Cl("++"))) # shorthand for metric signature (1, 1)
2-component Multivector{Cl("++"), 1, SVector{2, Int64}}:
 1 v1
 1 v2
```


### [The metric signature interface](@id sig-interface)


The metric signature type parameter may be any `isbits` value satisying the following interface.
As well as defining the geometric algebra, the signature is used to specify basis blade labels, the default array type for multivector components, and other metadata.

| Required methods | Description |
|:-----------------|:------------|
| `canonical_signature(sig)` | Canonical representation as a tuple (e.g., `Cl(1,3)` becomes `(1,-1,-1,-1)`).

| Optional methods | Description | Default
|:-----------------|:------------|:-------
| `dimension(sig)` | The dimension of the underlying vector space | `length(canonical_signature(sig))`
| `basis_vector_square(sig, i)` | The scalar square of the `i`th basis vector | `canonical_signature(sig)[i]`
| `show_signature(io, sig)` | Print the metric signature in a compact human-readable form | `show(io, sig)`
| `show_basis_blade(io, sig, indices)` | Print a basis blade with the given indices (e.g., `v12` or `𝒆₁∧𝒆₂`) | prints like `v12`
| `bits_to_indices(sig, bits)` | Define the order of indices for a basis blade (must also implement `basis_blade_parity(sig, bits)` consistently) | always increasing like `v123`
| `componentstype(sig, N)` | Preferred array type for `Multivector{sig}` components. (E.g., `Vector`, `MVector`, `SparseVector`, etc.) | `SVector` if small, `Vector` if large
| `use_symbolic_optim(sig)` | Whether to use symbolic code generation to optimise multivector products. | True for low dimensions


Below is an example of how one might define a “projectivised” signature which adds a projective dimension ``𝐯_0`` squaring to ``-1`` to any signature:
```@example ga
import GeometricAlgebra: canonical_signature, show_signature, show_basis_blade

struct ℙ{Sig} end
ℙ(sig) = ℙ{sig}()

canonical_signature(::ℙ{Sig}) where Sig = (-1, canonical_signature(Sig)...)
show_signature(io::IO, ::ℙ{Sig}) where Sig = print(io, "ℙ($Sig)")
show_basis_blade(io::IO, ::ℙ, indices::Vector) = print(io, "v", join(indices .- 1))

basis(ℙ(3))
```


## Symbolic Algebra and Code Generation

The package comes with a very minimal symbolic algebra module, `MiniCAS`, which is used for compile-time code simplification, but can also be used for (very simple) symbolic calculations.
For example, we can compute the product of two vectors symbolically as follows:

```jldoctest
julia> [Multivector{2,1}(:A), Multivector{2,2}(:B)]
2-element Vector{Multivector{2, K, Vector{GeometricAlgebra.MiniCAS.ProductNode{Expr}}} where K}:
 A[1]v1 + A[2]v2
 B[1]v12

julia> prod(ans)
2-component Multivector{2, 1, SVector{2, Any}}:
 -(A[2] * B[1]) v1
 A[1] * B[1]    v2
```

This makes it easy to optimize multivector operations: first perform the calculation symbolically and then compile the resulting analytic expression. By default, this compile-time optimization is enabled for most products (including the geometric, wedge and inner products in up to eight dimensions[^1]).
This optimisation is applied by prefixing method definitions with the internal [`@symbolic_optim`](@ref) macro.

[^1]: This can be changed on a per-algebra basis by defining methods for [`use_symbolic_optim()`](@ref).

Symbolic optimisation is also exposed through a user-facing macro [`@symbolicga`](@ref), inspired by the `@ga` macro in [serenity4/SymbolicGA.jl](https://github.com/serenity4/SymbolicGA.jl).
This is especially useful when you want use geometric algebra without manipulating `Multivector` types.

For example, using ``Cl(ℝ³)`` to represent homogeneous coordinates on the plane:
```jldoctest
julia> joinpoints(p, q) = @symbolicga 3 (p=1, q=1) p ∧ q Tuple
joinpoints (generic function with 1 method)

julia> intersectlines(p, q) = @symbolicga 3 (p=2, q=2) rdual(ldual(p) ∧ ldual(q)) Tuple
intersectlines (generic function with 1 method)

julia> L1 = joinpoints((1, 0, 1), (0, 1, 1)) # line y = 1 - x
(1, 1, -1)

julia> L2 = joinpoints((0, 0, 1), (1, 1, 1)) # line y = x
(0, -1, -1)

julia> intersectlines(L1, L2) # point (0.5, 0.5)
(-1, -1, -2)

```
The resulting methods are loop-free and allocation-free.
