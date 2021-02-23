using Documenter, GeometricAlgebra

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using GeometricAlgebra
end; recursive=true)

make() = makedocs(
	root=joinpath(dirname(pathof(GeometricAlgebra)), "..", "docs"),
	sitename="GeometricAlgebra.jl",
)

deploy() = deploydocs(
	repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

if "deploy" in ARGS
	make()
	deploy()
elseif "make" in ARGS
	make()
end