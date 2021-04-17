using GeometricAlgebra
using Documenter

DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using Revise, GeometricAlgebra
end; recursive=true)

doctest(GeometricAlgebra)