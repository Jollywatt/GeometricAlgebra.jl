using Documenter

DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using GeometricAlgebra
end; recursive=true)

doctest(GeometricAlgebra)