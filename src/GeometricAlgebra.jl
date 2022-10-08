"""
# GeometricAlgebra

Implements multivector (or ``k``-vector) types from geometric algebra (a.k.a. Clifford algebra).

## Exported Types

```
                   AbstractMultivector{Sig}
                     /                  \\
   HomogeneousMultivector{Sig,K}    Multivector{Sig,S}
       /               \\                             
Blade{Sig,K,T}   KVector{Sig,K,S}                
                                                   
                 ╰─── CompositeMultivector{Sig,S} ───╯
```

See [`basis`](@ref) and [`@basis`](@ref) to get started.

"""
module GeometricAlgebra

using StaticArrays, SparseArrays
using Combinatorics: permutations, powerset
using SymbolicUtils
using PrettyTables: pretty_table

export AbstractMultivector, HomogeneousMultivector, Blade, KVector, Multivector
export basis, @basis, @basisall, cayleytable
export signature, dimension, grade, scalar, ncomponents, blades
export geometric_prod, scalar_prod, ⊙, wedge, ∧, inner, ⋅, lcontract, ⨼, rcontract, ⨽, clifford_conj, var"'ᶜ"
export reversion, involution
export vector_repr, matrix_repr
export Cl

include("bits.jl")
include("types.jl")
include("utils.jl")
include("signatures.jl")
include("algebra.jl")
include("special.jl")
include("generated.jl")
include("show.jl")


end # module GeometricAlgebra
