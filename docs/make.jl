#= Run this script interactively: `julia -i make.jl`
... or with arguments `julia make.jl [fix|make|deploy]` =#

using Documenter, Revise, GeometricAlgebra

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using Revise, GeometricAlgebra
end; recursive=true)

make() = makedocs(
	root=joinpath(dirname(pathof(GeometricAlgebra)), "..", "docs"),
	sitename="GeometricAlgebra.jl",
)

deploy() = deploydocs(
	repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

fix() = begin
	Revise.revise()
	doctest(GeometricAlgebra, fix=true)
end

if isempty(ARGS)
	@info """Run this script in interactive mode, and call any of the following:
		`make()` - build documentation locally
		`deploy()` - build and deploy to github
		`fix()` - fix doctests
		Keep the session alive; changes will be revised and successive runs will be faster.
		"""
end

if "fix" in ARGS
	fix()
end

if "deploy" in ARGS
	make()
	deploy()
elseif "make" in ARGS
	make()
end
