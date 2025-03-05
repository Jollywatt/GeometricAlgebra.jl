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
         AbstractMultivector{Sig}
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

The metric signature type parameter `Sig` defines the dimension of the geometric algebra and the norms of its standard orthonormal basis vectors.
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
2-component Multivector{Cl("++"), 1, MVector{2, Int64}}:
 1 v1
 1 v2
```


### [The metric signature interface](@id sig-interface)


The metric signature type parameter may be any `isbits` value satisying the following interface.
As well as defining the geometric algebra, the signature is used to specify basis blade labels, the default array type for multivector components, and other metadata.

| Required methods | Description |
|:-----------------|:------------|
| `dimension(sig)` | The dimension of the underlying vector space, or number of basis vectors.
| `basis_vector_square(sig, i)` | The scalar square of the `i`th basis vector. |

| Optional methods | Description |
|:-----------------|:------------|
| `show_signature(io, sig)` | Show the metric signature in a compact human-readable form.
| `show_basis_blade(io, sig, indices)` | Print a basis blade with the given indices (e.g., `v12` or `𝒆₁∧𝒆₂`).
| `bits_to_indices(sig, bits)` | Define display order of indices for a basis blade (must also implement `basis_blade_parity(sig, bits)` consistently).
| `componentstype(sig, N)` | Preferred array type for `Multivector{sig}` components. (E.g., `Vector`, `MVector`, `SparseVector`, etc.)
| `use_symbolic_optim(sig)` | Whether to use symbolic code generation to optimise multivector products. (Default is true for low dimensions.)


Below is an example of how one might define a “projectivised” signature which adds a projective dimension ``𝐯_0`` squaring to ``-1`` to any signature:
```@example ga
import GeometricAlgebra: dimension, basis_vector_square, show_signature, show_basis_blade

struct ℙ{Sig} end
ℙ(sig) = ℙ{sig}()

dimension(::ℙ{Sig}) where Sig = dimension(Sig) + 1
basis_vector_square(::ℙ{Sig}, i) where Sig = i == 1 ? -1 : basis_vector_square(Sig, i - 1)
show_signature(io::IO, ::ℙ{Sig}) where Sig = print(io, "ℙ($Sig)")

show_basis_blade(io::IO, ::ℙ, indices::Vector) = print(io, "v", join(indices .- 1))

basis(ℙ(3)) |> sum
```


## Symbolic Algebra and Code Generation

Thanks to the wonderful [`SymbolicUtils`](https://symbolicutils.juliasymbolics.org/) package, the same code originally written for numerical multivectors readily works with symbolic components.
For example, we can compute the product of two vectors symbolically as follows:

```jldoctest
julia> GeometricAlgebra.make_symbolic.(Multivector{2,1}, [:A, :B])
2-element Vector{Multivector{2, 1, Vector{GeometricAlgebra.MiniCAS.ProductNode{GeometricAlgebra.MiniCAS.IndexNode{1}}}}}:
 A[1]v1 + A[2]v2
 B[1]v1 + B[2]v2

julia> prod(ans)
2-component Multivector{2, 0:2:2, MVector{2, GeometricAlgebra.MiniCAS.SumNode{GeometricAlgebra.MiniCAS.IndexNode{1}, Int64}}}:
 B[1] * A[1] + B[2] * A[2]
 -(B[1] * A[2]) + B[2] * A[1] v12
```

This makes it easy to optimize multivector operations: first perform the calculation symbolically and then compile the resulting analytic expression. By default, this optimization is enabled for most products (including the geometric, wedge and inner products in up to eight dimensions[^1]).
This is done by prefixing method definitions with the internal [`@symbolic_optim`](@ref) macro.

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
