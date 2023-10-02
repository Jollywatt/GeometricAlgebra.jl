using GeometricAlgebra:
	SingletonVector,
	superscript,
	subscript

@testset "SingletonVector" begin
	@test collect(SingletonVector(42, 2, 5)) == [0, 42, 0, 0, 0]
	@test size(SingletonVector(42, 5, 10_000)) == (10_000,)

	sv = SingletonVector(1, 5, 10_000)
	@test 7*sv == SingletonVector(7, 5, 10_000)
	@test sv/4 == SingletonVector(0.25, 5, 10_000)
	@test !iszero(sv)
	@test iszero(0sv)

	nan = SingletonVector("nan", 2, 5)
	@test nan[1] === 0 && nan[2] === "nan"
	@test eltype(nan) === Any
	@test collect(nan) == [0, "nan", 0, 0, 0]
end

@testset "super/subscript" begin
	@test superscript(42) == "⁴²"
	@test subscript(-808) == "₋₈₀₈"
	@test collect(subscript(10123456789))[2:end] == GeometricAlgebra.SUBSCRIPT_DIGITS
end