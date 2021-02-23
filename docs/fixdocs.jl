#= should be run from `./GeometricAlgebra.jl/` interactively =#

using Documenter
using Revise, GeometricAlgebra

DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using Revise, GeometricAlgebra
end; recursive=true)

main() = while true
	print("Press return to fix docs: ")
	isempty(readline(stdin)) || break
	Revise.revise()
	doctest(GeometricAlgebra, fix=true)
end

main()