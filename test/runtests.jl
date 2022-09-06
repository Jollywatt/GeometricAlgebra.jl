#= Run this script interactively: `julia --project -i runtests.jl`
... or with arguments `julia --project runtests.jl [testfiles...]` =#

cd(dirname(PROGRAM_FILE))
using Pkg; Pkg.activate(".")

using Random, Test, Coverage
using Revise, Multivectors


alltests() = setdiff(filter(endswith(".jl"), readdir()), [basename(PROGRAM_FILE)])

function test(files=alltests())
	
	@testset "$file" for file in files
		include(file)
	end
	nothing
end

function coverage()
	run(`julia --code-coverage --project runtests.jl`)
	score()
end

function score()
	local percentage
	with_logger(NullLogger()) do
		coverage = process_folder(joinpath(project_root, "src"))
		covered, total = get_summary(coverage)
		percentage = round(100covered/total, digits=2)
	end
	@info "Coverage: $percentage%"
end

clear_coverage_files() = run(`bash -c "rm $project_root/**/*.cov"`)


if isempty(ARGS)
	@info """Run this script in interactive mode, and call:
		 - `test("testfile1", "testfile2", ...)`
		   (without `.jl` extensions) to run tests
		 - `test()` to run all tests
		 - `coverage()` to analyse code coverage
		 - `score()` to show coverage score from last run
		 - `clear_coverage_files()` to remove `*.cov` files
		Keep the session alive; changes will be revised and successive runs will be faster.
		"""
	!isinteractive() && test()
elseif "--code-coverage" in ARGS
	coverage()
else
	test(ARGS...)
end
