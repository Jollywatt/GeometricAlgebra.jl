#= Run this script interactively: `julia --project=docs/ -i make.jl`
... or with arguments `julia --project=docs/ make.jl [test|fix|make|deploy]` =#

using Documenter, Revise, GeometricAlgebra

project_root = dirname(dirname(pathof(GeometricAlgebra)))

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using Revise, GeometricAlgebra
end; recursive=true)

make() = makedocs(
	root=joinpath(project_root, "docs"),
	sitename="GeometricAlgebra.jl",
	pages=[
		"Home"=>"index.md",
	],
)

deploy() = deploydocs(
	repo = "github.com/Jollywatt/GeometricAlgebra.jl.git",
)

test() = doctest(GeometricAlgebra)

fix() = begin
	Revise.revise()
	doctest(GeometricAlgebra, fix=true)
end

if isempty(ARGS)
	@info """Run this script in interactive mode, and call:
		 - `make()` to build documentation locally
		 - `deploy()` to build and deploy to github
		 - `test()` to run doctests
		 - `fix()` to fix doctests
		Keep the session alive; changes will be revised and successive runs will be faster.
		Alternatively, run this script passing a command (without parentheses) as an argument.
		"""
end

if "test" in ARGS
	test()
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
