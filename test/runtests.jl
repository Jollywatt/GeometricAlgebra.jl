#= Run this script interactively: `julia -i runtests.jl`
... or with arguments `julia runtests.jl [testfiles...]` =#

using Random, Test, Coverage, Documenter
using Revise, GeometricAlgebra

const project_root = pathof(GeometricAlgebra) |> dirname |> dirname

# apply setup code to all doctests in doc strings
DocMeta.setdocmeta!(GeometricAlgebra, :DocTestSetup, quote
	using Revise
	using GeometricAlgebra
	using GeometricAlgebra.MiniCAS
end; recursive=true)

function alltests(dir=dirname(@__FILE__))
	testfiles = String[]
	for (root, dirs, files) in walkdir(dir)
		f = filter(files) do file
			endswith(file, ".jl") && file != @__FILE__
		end
		append!(testfiles, relpath.(joinpath.(root, f), pwd()))
		for dir in dirs
			append!(testfiles, alltests(dir))
		end
	end
	testfiles
end

test(files::String...) = test(files)
function test(files=alltests())
	@testset "$file" for file in files
		include(file)
	end
	nothing
end

function coverage()
	clear_coverage_files()
	run(`julia --code-coverage --project runtests.jl`)
	score()
end

function score()
	coverage = process_folder(joinpath(project_root, "src"))
	covered, total = get_summary(coverage)
	percentage = round(100covered/total, digits=2)
	@info "Overall coverage: $percentage%"
end

clear_coverage_files() = run(`bash -c "rm -f $project_root/**/*.cov"`)


if isempty(ARGS)
	@info """Run this script in interactive mode, and call:
		 - `test("testfile1.jl", "testfile2.jl", ...)`
		 - `test()` to run all tests
		 - `coverage()` to analyse code coverage
		 - `score()` to show coverage score from last run
		 - `clear_coverage_files()` to remove `*.cov` files
		Keep the session alive; changes will be revised and successive runs will be faster.
		"""
	if !isinteractive()
		test()
		VERSION >= v"1.11" && doctest(GeometricAlgebra) # doctests are sensitive to Julia version
	end
elseif "--code-coverage" in ARGS
	coverage()
else
	test(ARGS...)
end
