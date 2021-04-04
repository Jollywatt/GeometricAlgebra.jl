#= Run this script interactively: `julia --project=test/ -i runtests.jl`
... or with arguments `julia --project=test/ runtests.jl` =#

using Test, TestSetExtensions, Random
using Revise, GeometricAlgebra

function runtests(files...)
	@eval @testset ExtendedTestSet "GeometricAlgebra" begin
		@includetests $files
	end
	nothing
end

cd(joinpath(dirname(pathof(GeometricAlgebra)), "..", "test"))

if isempty(ARGS)
	@info """Run this script in interactive mode, and call
		`runtests("testfile1", "testfile2", ...)`
		(without `.jl` extensions) to run tests.
		Call `runtests()` to run all tests.
		Keep the session alive; changes will be revised and successive runs will be faster.
		"""
	!isinteractive() && runtests()
else
	runtests(ARGS...)
end
