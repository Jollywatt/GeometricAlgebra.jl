using Documenter, GeometricAlgebra

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using GeometricAlgebra
end; recursive=true)

cd(pathof(GeometricAlgebra))

make() = makedocs(
	root="docs",
	sitename="GeometricAlgebra.jl",
)

deploy() = deploydocs(
	repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

if "deploy" in ARGS
	make()
	deploy()
end
