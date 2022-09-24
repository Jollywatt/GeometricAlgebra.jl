var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [GeometricAlgebra]","category":"page"},{"location":"reference/#GeometricAlgebra.GeometricAlgebra","page":"Reference","title":"GeometricAlgebra.GeometricAlgebra","text":"GeometricAlgebra\n\nImplements multivector (or k-vector) types from geometric algebra (a.k.a. Clifford algebra).\n\nExported Types\n\n                   AbstractMultivector{Sig}\n                     /                  \\\n   HomogeneousMultivector{Sig,K}    MixedMultivector{Sig,S}\n       /                \\                             \nBlade{Sig,K,T}    Multivector{Sig,K,S}                \n                                                   \n                  ╰───── CompositeMultivector{Sig,S} ─────╯\n\nSee basis and @basis to get started.\n\n\n\n\n\n","category":"module"},{"location":"reference/#GeometricAlgebra.AbstractMultivector","page":"Reference","title":"GeometricAlgebra.AbstractMultivector","text":"AbstractMultivector{Sig}\n\nSupertype of all elements in the geometric algebra defined by the metric signature Sig.\n\nSubtypes\n\n                   AbstractMultivector{Sig}\n                     /                  \\\n   HomogeneousMultivector{Sig,K}    MixedMultivector{Sig,S}\n       /                \\                             \nBlade{Sig,K,T}    Multivector{Sig,K,S}                \n                                                   \n                  ╰───── CompositeMultivector{Sig,S} ─────╯\n\nBlade: a scalar multiple of a wedge product of orthogonal basis vectors.\nMultivector: a homogeneous multivector; a sum of same-grade blades.\nMixedMultivector: an inhomogeneous multivector. All elements in a geometric  algebra can be represented as this type (though not most efficiently).\n\nnote: Note\nThe mathematical definition of a k-blade is the wedge product of k vectors, not necessarily basis vectors. Thus, not all k-blades are representable as a Blade, but are always representable as a sum of Blades, or as a Multivector.\n\nType Parameters\n\nSig: The metric signature which defines the geometric algebra. This can be any  all-bits value which satisfies the metric signature interface.  For example, (1, 1, 1) or EuclideanSignature(3) both  define the standard geometric algebra of ℝ^3.\nT: The numerical type of the coefficient of a Blade.\nK: An Int specifying the grade of a HomogeneousMultivector.\nS: The storage type of the components of a CompositeMultivector. This is  assumed to be mutable, and is usually a subtype of Vector, MVector or SparseVector.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.BitPermutations","page":"Reference","title":"GeometricAlgebra.BitPermutations","text":"BitPermutations{T}(n)\n\nInfinite iterator returning all unsigned integers of type T, in ascending order, for which Base.count_ones is n.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.Blade","page":"Reference","title":"GeometricAlgebra.Blade","text":"Blade{Sig,K,T} <: HomogeneousMultivector{Sig,K}\n\nA blade of grade K with basis blade bits and scalar coefficient of type T.\n\nParameters\n\nSig: Metric signature defining the geometric algebra, retrieved with signature().\nK: Grade of the blade, equal to count_ones(bits), retrieved with grade().\nT: Numerical type of the scalar coefficient, retrieved with eltype().\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.Blade-Union{Tuple{Pair}, Tuple{Sig}} where Sig","page":"Reference","title":"GeometricAlgebra.Blade","text":"Blade{Sig}(bits, coeff)\nBlade{Sig}(bits => coeff)\n\nBasis blade with indices encoded by bits and scalar coefficient coeff.\n\nExamples\n\njulia> Blade{3}(0b110 => 42) # a grade 2 blade in 3 dimensions\nBlade{3, 2, Int64}:\n 42 v23\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.Cl","page":"Reference","title":"GeometricAlgebra.Cl","text":"Cl(p, q=0, r=0)\n\nMetric signature where p, q and r is the number of basis vectors of norm +1, -1 and 0, respectively.\n\nExamples\n\njulia> basis(Cl(1,3))\n4-element Vector{Blade{Cl(1,3), 1, Int64}}:\n 1v1\n 1v2\n 1v3\n 1v4\n\njulia> ans .^ 2\n4-element Vector{Blade{Cl(1,3), 0, Int64}}:\n  1\n -1\n -1\n -1\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.HomogeneousMultivector","page":"Reference","title":"GeometricAlgebra.HomogeneousMultivector","text":"HomogeneousMultivector{Sig,K} <: AbstractMultivector{Sig}\n\nAbstract supertype of Blade and Multivector.\n\nParameters\n\nSig: Metric signature defining the geometric algebra, retrieved with signature().\nK: Grade of the blade or multivector, retrieved with grade().\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.MixedMultivector","page":"Reference","title":"GeometricAlgebra.MixedMultivector","text":"MixedMultivector{Sig,S} <: AbstractMultivector{Sig}\n\nA (possibly inhomogeneous) multivector.\n\nAll elements of a geometric algebra are representable as a MixedMultivector.\n\nParameters\n\nSig: Metric signature defining the geometric algebra, retrieved with signature().\nS: Storage type of the multivector components, usually a subtype of AbstractVector.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.MixedMultivector-Union{Tuple{S}, Tuple{Sig}} where {Sig, S}","page":"Reference","title":"GeometricAlgebra.MixedMultivector","text":"MixedMultivector{Sig}(comps)\n\nInhomogeneous multivector with components vector comps. The components are ordered first by grade then lexicographically (see GeometricAlgebra.mmv_bits).\n\nExamples\n\njulia> MixedMultivector{3}(1:2^3)\n8-component MixedMultivector{3, UnitRange{Int64}}:\n 1\n 2 v1 + 3 v2 + 4 v3\n 5 v12 + 6 v13 + 7 v23\n 8 v123\n\njulia> grade(ans, 1)\n3-component Multivector{3, 1, UnitRange{Int64}}:\n 2 v1\n 3 v2\n 4 v3\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.Multivector","page":"Reference","title":"GeometricAlgebra.Multivector","text":"Multivector{Sig,K,S} <: HomogeneousMultivector{Sig,K}\n\nA homogeneous multivector of grade K with storage type S.\n\nParameters\n\nSig: Metric signature defining the geometric algebra, retrieved with signature().\nK: Grade of the multivector, retrieved with grade().\nS: Storage type of the multivector components, usually a subtype of AbstractVector.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GeometricAlgebra.Multivector-Union{Tuple{S}, Tuple{K}, Tuple{Sig}} where {Sig, K, S}","page":"Reference","title":"GeometricAlgebra.Multivector","text":"Multivector{Sig,K}(comps)\n\nMultivector of grade K with components vector comps.\n\nExamples\n\njulia> Multivector{3,2}(1:3) # 3D bivector\n3-component Multivector{3, 2, UnitRange{Int64}}:\n 1 v12\n 2 v13\n 3 v23\n\n\n\n\n\n","category":"method"},{"location":"reference/#Base.eltype-Union{Tuple{Union{Type{var\"#s4\"}, var\"#s4\"} where var\"#s4\"<:(Blade{Sig, K, T} where {Sig, K})}, Tuple{T}} where T","page":"Reference","title":"Base.eltype","text":"eltype(::AbstractMultivector)\n\nThe numerical type of the components of a multivector instance or type.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.basis-Tuple{Any}","page":"Reference","title":"GeometricAlgebra.basis","text":"basis(sig; grade=1)\n\nReturn basis blades for the geometric algebra defined by the metric signature sig.\n\nExamples\n\njulia> basis(3)\n3-element Vector{Blade{3, 1, Int64}}:\n 1v1\n 1v2\n 1v3\n\njulia> basis(Cl(1,3); grade=2)\n6-element Vector{Blade{Cl(1,3), 2, Int64}}:\n 1v12\n 1v13\n 1v23\n 1v14\n 1v24\n 1v34\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.bits_of_grade-Tuple{Any}","page":"Reference","title":"GeometricAlgebra.bits_of_grade","text":"bits_of_grade(k[, dim])\n\nGenerate basis blade bits of grade k in ascending order. Yields all basis blades in the dimension dim, if given, otherwise iterate indefinitely.\n\nExamples\n\njulia> GeometricAlgebra.bits_of_grade(2, 4) .|> UInt8 .|> bitstring\n6-element Vector{String}:\n \"00000011\"\n \"00000101\"\n \"00000110\"\n \"00001001\"\n \"00001010\"\n \"00001100\"\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.bits_to_indices-Tuple{Unsigned}","page":"Reference","title":"GeometricAlgebra.bits_to_indices","text":"bits_to_indices(bits)\n\nReturn the positions of the ones in the unsigned integer bits.\n\nUsed to convert between representations of a unit blade. Inverse of indices_to_bits.\n\nExamples\n\njulia> GeometricAlgebra.bits_to_indices(0b1001101)\n4-element Vector{Int64}:\n 1\n 3\n 4\n 7\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.bits_to_mmv_index-Tuple{Unsigned, Any}","page":"Reference","title":"GeometricAlgebra.bits_to_mmv_index","text":"bits_to_mmv_index(bits::Unsigned)\n\nConvert a unit blade bits to a linear index for accessing components of a MixedMultivector. \n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.bits_to_mv_index-Tuple{Unsigned}","page":"Reference","title":"GeometricAlgebra.bits_to_mv_index","text":"bits_to_mv_index(bits::Unsigned)\n\nConvert a unit blade bits to a linear index for accessing components of a Multivector. \n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.componentstype-Tuple{Any, Any, Any}","page":"Reference","title":"GeometricAlgebra.componentstype","text":"componentstype(sig, N, T)\n\nArray type to use the components for multivectors of signature sig. The resulting type should be able to store N components (in the case of a fixed-size array) of element type T.\n\nThis is used when converting a Blade into a CompositeMultivector to determine the type of the components array.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.dimension-Union{Tuple{AbstractMultivector{Sig}}, Tuple{Sig}} where Sig","page":"Reference","title":"GeometricAlgebra.dimension","text":"dimension(::AbstractMultivector)\n\nThe dimension of the underlying vector space of the geometric algebra. See ncomponents for the dimension of the algebra (i.e., the number of independent components of a general multivector).\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.factor_from_squares-Tuple{Any, Unsigned}","page":"Reference","title":"GeometricAlgebra.factor_from_squares","text":"factor_from_squares(sig, bits::Unsigned)\n\nCompute the overall factor arising from the geometric product between repeated basis vectors.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.generated_multivector_function-Tuple{Any, Vararg{Any}}","page":"Reference","title":"GeometricAlgebra.generated_multivector_function","text":"generated_multivector_function(f, a, b, ...)\n\nTrace evaluation of f(a, b, ...)::CompositeMultivector on symbolic versions of each AbstractMultivector instance or type a, b, ..., returning an expression body suitable for a @generated function.\n\nThe names of the arguments to be passed to the function body are the literal symbols :a, :b, etc.\n\nExamples\n\njulia> u, v = Multivector.(basis(2))\n2-element Vector{Multivector{2, 1, Vector{Int64}}}:\n 1v1\n 1v2\n\njulia> using MacroTools: prettify\n\njulia> ex = GeometricAlgebra.generated_multivector_function(*, u, v) |> prettify\n:(let a = a.components, b = b.components\n      (MixedMultivector{2})([a[1] * b[1] + a[2] * b[2], 0, 0, a[1] * b[2] + (-1 * a[2]) * b[1]])\n  end)\n\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.geometric_prod_bits-Tuple{Any, Unsigned, Unsigned}","page":"Reference","title":"GeometricAlgebra.geometric_prod_bits","text":"geometric_prod_bits(sig, a::Unsigned, b::Unsigned)\n\nCompute the geometric product between unit blades. Returns a tuple of the overall scalar factor and the resulting unit blade.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.grade-Union{Tuple{Union{Type{var\"#s4\"}, var\"#s4\"} where var\"#s4\"<:HomogeneousMultivector{Sig, K}}, Tuple{K}, Tuple{Sig}} where {Sig, K}","page":"Reference","title":"GeometricAlgebra.grade","text":"grade(::HomogeneousMultivector{Sig,K}) -> K\n\nThe grade of a homogeneous multivector (a Blade or Multivector) instance or type.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.indices_to_bits-Tuple{Any}","page":"Reference","title":"GeometricAlgebra.indices_to_bits","text":"indices_to_bits(indices)\n\nCreate unsigned integer with bits at the positions given in the vector indices.\n\nUsed to convert between representations of a unit blade. Inverse of bits_to_indices.\n\nwarning: Warning\nProduces incorrect results if elements of indices are greater than the number of bits in bits_scalar() <: Unsigned.\n\nExamples\n\njulia> GeometricAlgebra.indices_to_bits([1, 2, 5]) |> UInt16 |> bitstring\n\"0000000000010011\"\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.mmv_bits-Union{Tuple{Val{N}}, Tuple{N}} where N","page":"Reference","title":"GeometricAlgebra.mmv_bits","text":"mmv_bits(::Val{N})\n\nVector of unit blades corresponding to components of an N-dimensional MixedMultivector. Mixed multivector components are ordered first by grade then by numerical value of the unit blade, so that the grade K components are contiguous and given by mv_bits(Val(N), Val(K))\n\nExamples\n\njulia> GeometricAlgebra.mmv_bits(Val(3)) .|> UInt8 .|> bitstring\n8-element Vector{String}:\n \"00000000\"\n \"00000001\"\n \"00000010\"\n \"00000100\"\n \"00000011\"\n \"00000101\"\n \"00000110\"\n \"00000111\"\n\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.mv_bits-Union{Tuple{K}, Tuple{N}, Tuple{Val{N}, Val{K}}} where {N, K}","page":"Reference","title":"GeometricAlgebra.mv_bits","text":"mv_bits(::Val{N}, ::Val{K})\n\nVector of unit blades corresponding to components of an N-dimensional Multivector of grade K. Multivector components are sorted by the numerical value of the unit blade.\n\nExamples\n\njulia> GeometricAlgebra.mv_bits(Val(4), Val(2)) .|> UInt8 .|> bitstring\n6-element Vector{String}:\n \"00000011\"\n \"00000101\"\n \"00000110\"\n \"00001001\"\n \"00001010\"\n \"00001100\"\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.ncomponents-Tuple{Any}","page":"Reference","title":"GeometricAlgebra.ncomponents","text":"ncomponents(sig)\nncomponents(sig, k)\n\nDimension of (the grade-k subspace of) the geometric algebra of metric signature sig, viewed as a vector space.\n\nIf the dimension of the underlying vector space in n, then the algebra is 2^n-dimensional, and its grade-k subspace binomnk-dimensional.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.ncomponents-Union{Tuple{Union{Type{var\"#s4\"}, var\"#s4\"} where var\"#s4\"<:(Multivector{Sig, K})}, Tuple{K}, Tuple{Sig}} where {Sig, K}","page":"Reference","title":"GeometricAlgebra.ncomponents","text":"ncomponents(::CompositeMultivector)\n\nNumber of independent components of a multivector instance or type.\n\nIn n dimensions, this is binomnk for a Multivector and 2^n for a MixedMultivector.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.next_bit_permutation-Tuple{Unsigned}","page":"Reference","title":"GeometricAlgebra.next_bit_permutation","text":"Return the smallest uint larger than the one given which has the same number of binary ones. Algorithm is Gosper’s hack.\n\nExamples\n\njulia> GeometricAlgebra.next_bit_permutation(0b1011) |> bitstring\n\"00001101\"\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.show_basis_blade-Tuple{Any, Any, Any}","page":"Reference","title":"GeometricAlgebra.show_basis_blade","text":"show_basis_blade(io, sig, indices)\n\nShow the basis blade which is the product of the unit vectors in indices in a geometric algebra defined by sig. Methods dispatching on sig should be added to customise basis blade labels for particular algebras.\n\nThe fallback method is:\n\nshow_basis_blade(io, sig, indices) = printstyled(io, \"v\"*join(string.(indices)); bold=true)\n\nExamples\n\njulia> GeometricAlgebra.show_basis_blade(stdout, (1, 1, 1), [1, 3])\nv13\n\njulia> using GeometricAlgebra: subscript\n\njulia> GeometricAlgebra.show_basis_blade(io, sig, indices) = print(io, join(\"𝒆\".*subscript.(indices), \"∧\"))\n\njulia> prod(basis(4))\nBlade{⟨++++⟩, 4, Int64} of grade 4:\n 1 𝒆₁∧𝒆₂∧𝒆₃∧𝒆₄\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.show_blade-Tuple{IO, Blade}","page":"Reference","title":"GeometricAlgebra.show_blade","text":"Display blade with parentheses surrounding coefficient if necessary.\n\njulia> GeometricAlgebra.show_blade(stdout, Blade{(x=1,)}(0b1 => 1 + im))\n(1+1im) x\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.show_mixedmultivector-Tuple{IO, MixedMultivector}","page":"Reference","title":"GeometricAlgebra.show_mixedmultivector","text":"Display an inhomogeneous MixedMultivector with each grade on a new line.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.show_multivector-Tuple{IO, Multivector}","page":"Reference","title":"GeometricAlgebra.show_multivector","text":"Display homogeneous multivector components as a column of blades, with coefficients and blades aligned using the native alignment mechanism.\n\njulia> a = Multivector{(1,1,1),1}([1e3, 1, 1e-3]);\n\njulia> GeometricAlgebra.show_multivector(stdout, a)\n1000.0   v1\n   1.0   v2\n   0.001 v3\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.show_signature-Tuple{Any, Any}","page":"Reference","title":"GeometricAlgebra.show_signature","text":"show_signature(io, sig)\n\nPretty-print the metric signature sig.\n\nThis is used to display the metric signature type parameter in AbstractMultivector subtypes to reduce visual noise. Methods may optionally be added for user-defined metric signatures, in a similar fashion to Base.show.\n\nExamples\n\njulia> sig = (+1,-1,-1,-1)\n(1, -1, -1, -1)\n\njulia> GeometricAlgebra.show_signature(stdout, sig)\n⟨+---⟩\n\njulia> Blade{sig}\nBlade{⟨+---⟩} (pretty-printed Blade{(1, -1, -1, -1)})\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.sign_from_swaps-Tuple{Unsigned, Unsigned}","page":"Reference","title":"GeometricAlgebra.sign_from_swaps","text":"Compute sign flips of blade product due to transposing basis vectors into sorted order. (The full sign of the product will also depend on the basis norms.)\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.signature-Union{Tuple{Union{Type{var\"#s4\"}, var\"#s4\"} where var\"#s4\"<:AbstractMultivector{Sig}}, Tuple{Sig}} where Sig","page":"Reference","title":"GeometricAlgebra.signature","text":"signature(::AbstractMultivector{Sig}) -> Sig\n\nThe metric signature type parameter of the multivector instance or type.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GeometricAlgebra.@basis-Tuple{Any}","page":"Reference","title":"GeometricAlgebra.@basis","text":"@basis sig\n\nPopulate namespace with basis blades of every grade in the geometric algebra with metric signature sig.\n\nSee also @basisall.\n\nExamples\n\njulia> @basis 3\n[ Info: Defined basis blades v, v1, v2, v3, v12, v13, v23, v123\n\njulia> 1v2 + 3v12\n8-component MixedMultivector{3, Vector{Int64}}:\n 1 v2\n 3 v12\n\n\n\n\n\n","category":"macro"},{"location":"reference/#GeometricAlgebra.@basisall-Tuple{Any}","page":"Reference","title":"GeometricAlgebra.@basisall","text":"@basisall sig\n\nSimilarly to @basis, populate namespace with basis blades, but include all permutations of each blade.\n\nwarning: Warning\nThis defines 2^n variables for an n dimensional signature!\n\nExamples\n\njulia> @basisall (+1,-1)\n[ Info: Defined basis blades v, v1, v2, v12, v21\n\njulia> v12 == -v21\ntrue\n\n\n\n\n\n","category":"macro"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"CurrentModule = GeometricAlgebra\nDocTestSetup = quote\n\tusing GeometricAlgebra\nend","category":"page"},{"location":"#GeometricAlgebra","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"","category":"section"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"GeometricAlgebra.jl implements basic types for working with geometric (or Clifford) algebras.","category":"page"},{"location":"#Quick-Start","page":"GeometricAlgebra","title":"Quick Start","text":"","category":"section"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Construct multivectors by providing the metric signature and grade as type parameters:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"julia> u = Multivector{(1,1,1),1}([1, -1, 0])\n3-component Multivector{⟨+++⟩, 1, Vector{Int64}}:\n  1 v1\n -1 v2\n  0 v3\n\njulia> v = Multivector{(1,1,1),2}(1:3)\n3-component Multivector{⟨+++⟩, 2, UnitRange{Int64}}:\n 1 v12\n 2 v13\n 3 v23\n\njulia> u*v + π\n8-component MixedMultivector{⟨+++⟩, Vector{Float64}}:\n 3.14159\n 1.0 v1 + 1.0 v2 + -1.0 v3\n 5.0 v123","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"You may also obtain an orthonormal basis for a metric signature:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"julia> v = basis((-1,+1,+1,+1))\n4-element Vector{Blade{⟨-+++⟩, 1, Int64}}:\n 1v1\n 1v2\n 1v3\n 1v4\n\njulia> exp(10000*2π*v[3]v[4])\n16-component MixedMultivector{⟨-+++⟩, Vector{Float64}}:\n 1.0\n -9.71365e-13 v34","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Macros are provided for interactive use:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"julia> @basis \"+---\"\n[ Info: Defined basis blades v, v1, v2, v3, v4, v12, v13, v14, v23, v24, v34, v123, v124, v134, v234, v1234\n\njulia> v1234^2\nBlade{⟨+---⟩, 0, Int64}:\n -1","category":"page"},{"location":"#Multivector-Types","page":"GeometricAlgebra","title":"Multivector Types","text":"","category":"section"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"There are three concrete types for representing elements in a geometric algebra, arranged in the following type hierarchy:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"                   AbstractMultivector{Sig}\n                     /                  \\\n   HomogeneousMultivector{Sig,K}    MixedMultivector{Sig,S}\n       /                \\                             \nBlade{Sig,K,T}    Multivector{Sig,K,S}                \n                                                   \n                  ╰───── CompositeMultivector{Sig,S} ─────╯","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Blade: a scalar multiple of a wedge product of orthogonal basis vectors.\nMultivector: a homogeneous multivector; a sum of same-grade blades.\nMixedMultivector: an inhomogeneous multivector. All elements in a geometric  algebra can be represented as this type (though not most efficiently).","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"note: Note\nThe mathematical definition of a k-blade is the wedge product of k vectors, not necessarily basis vectors. Thus, not all k-blades are representable as a Blade, but are always representable as a sum of Blades, or a Multivector.","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"These types have up to three of type parameters:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Sig: The metric signature which defines the geometric algebra. This can be any  all-bits value which satisfies the metric signature interface.\nT: The numerical type of the coefficient of a Blade.\nK: An Int specifying the grade of a HomogeneousMultivector.\nS: The storage type of the components of a CompositeMultivector. This is  assumed to be mutable, and is usually a subtype of Vector, MVector or SparseVector.","category":"page"},{"location":"#Metric-Signatures","page":"GeometricAlgebra","title":"Metric Signatures","text":"","category":"section"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"The metric signature Sig defines the dimension of the geometric algebra and the norms of its standard orthonormal basis vectors.","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"julia> sig = (-1,+1,+1,+1)\n(-1, 1, 1, 1)\n\njulia> Multivector{sig,1}(1:4)\n4-component Multivector{⟨-+++⟩, 1, UnitRange{Int64}}:\n 1 v1\n 2 v2\n 3 v3\n 4 v4\n\njulia> signature(ans)\n(-1, 1, 1, 1)\n","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Additionally, the Sig type parameter carries metadata via multiple dispatch. This allows various behaviours to be customised for each signature, including the basis blade labels, the default array type of multivector components, and so on.","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Below is an example of how one might define a geometric algebra with specific behaviours:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"struct DiracGamma end\n\n# define the algebra\nGeometricAlgebra.dimension(::DiracGamma) = 4\nGeometricAlgebra.basis_vector_norm(::DiracGamma, i) = i > 1 ? -1 : +1\n\n# set the preferred component storage type\nusing StaticArrays\nGeometricAlgebra.componentstype(::DiracGamma, N, T) = MVector{N,T}\n\n# custom labels\nGeometricAlgebra.show_basis_blade(io, ::DiracGamma, indices) = print(io, join(\"γ\".*GeometricAlgebra.superscript.(indices)))\n\nbasis(DiracGamma())\n# output\n4-element Vector{Blade{DiracGamma(), 1, Int64}}:\n 1γ¹\n 1γ²\n 1γ³\n 1γ⁴","category":"page"},{"location":"#sig_interface","page":"GeometricAlgebra","title":"Metric signature interface","text":"","category":"section"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"The metric signature type parameter may be any isbits value satisying the following interface.","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Required methods Description\ndimension(sig) The dimension of the underlying vector space, or number of basis vectors.\nbasis_vector_norm(sig, i) The norm of the ith basis vector.","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Optional methods Description\nshow_signature(io, sig) Show the metric signature in a compact form.\nshow_basis_blade(io, sig, indices) Print the basis blade with the given indices.","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"The methods above are predefined for:","category":"page"},{"location":"","page":"GeometricAlgebra","title":"GeometricAlgebra","text":"Ints, defining a Euclidean metric of that dimension,\nTuples, defining the norms of each basis vector,\nNamedTuples, defining the norms and labels.","category":"page"}]
}