@testset "basis, @basis, @basisall" begin
	@test length(basis(10)) == 10
	@test grade(prod(basis(5))) == 5

	@test basis((1,1,1)) == basis("+++")

	let
		@basis "-+++"
		@test v1^2 == v234^2 == -1
	end

	let
		@basisall "-+++"
		@test v123 == v231 == v312
	end

	let
		@basisall (x=1, t=-1)
		@test x^2 == xt^2 == tx^2 == -t^2 == 1
	end

end

@testset "cayleytable" begin
	@test sprint(cayleytable, 1:5, +) == """
	 (↓) + (→) │ 1  2  3  4   5
	───────────┼────────────────
	         1 │ 2  3  4  5   6
	         2 │ 3  4  5  6   7
	         3 │ 4  5  6  7   8
	         4 │ 5  6  7  8   9
	         5 │ 6  7  8  9  10"""

	@test sprint(cayleytable, "+++", ∧) == """
	 (↓) ∧ (→) │    1 │   v1     v2    v3 │  v12    v13   v23 │ v123
	───────────┼──────┼───────────────────┼───────────────────┼──────
	         1 │    1 │   v1     v2    v3 │  v12    v13   v23 │ v123
	───────────┼──────┼───────────────────┼───────────────────┼──────
	        v1 │   v1 │    0    v12   v13 │    0      0  v123 │    0
	        v2 │   v2 │ -v12      0   v23 │    0  -v123     0 │    0
	        v3 │   v3 │ -v13   -v23     0 │ v123      0     0 │    0
	───────────┼──────┼───────────────────┼───────────────────┼──────
	       v12 │  v12 │    0      0  v123 │    0      0     0 │    0
	       v13 │  v13 │    0  -v123     0 │    0      0     0 │    0
	       v23 │  v23 │ v123      0     0 │    0      0     0 │    0
	───────────┼──────┼───────────────────┼───────────────────┼──────
	      v123 │ v123 │    0      0     0 │    0      0     0 │    0"""
end