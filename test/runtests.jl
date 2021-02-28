using Test, TestSetExtensions
using Revise, GeometricAlgebra

function runtests(files...)
	@eval @testset ExtendedTestSet "GeometricAlgebra" begin
		@includetests $files
	end
	nothing
end

if isempty(ARGS)
	@info """Run this script in intreactive mode, and call
		`@runtests "testfile1" "testfile2" ...`
		(without `.jl` extensions) to run tests.
		Keep the session alive; changes will be revised and successive runs will be faster.
		"""
else
	runtests(ARGS...)
end

