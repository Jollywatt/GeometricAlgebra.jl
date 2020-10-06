using Documenter, GeometricAlgebra

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using GeometricAlgebra
end; recursive=true)

@show readdir(dirname(pathof(GeometricAlgebra)))

make() = makedocs(
	root=joinpath(pathof(GeometricAlgebra), "docs"),
	sitename="GeometricAlgebra.jl",
)

deploy() = deploydocs(
	repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

if "deploy" in ARGS
	make()
	deploy()
end
