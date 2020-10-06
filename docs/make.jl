using Pkg, Documenter
using GeometricAlgebra

# run all doctests in doc strings afer using GeometricAlgebra
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using GeometricAlgebra
end; recursive=true)

function make(; kwargs...)
	makedocs(
		root=Pkg.dir("GeometricAlgebra", "docs"),
		sitename="GeometricAlgebra.jl",
		modules=[GeometricAlgebra];
		kwargs...)

	deploydocs(
		repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
		versions = ["stable" => "v^", "v#.#", devurl => devurl],
	)
end

make()
