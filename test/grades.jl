using GeometricAlgebra:
	unify_grades,
	resulting_multivector_type

@testset "unify_grades" begin
	# canonicalization
	@test unify_grades(3, 0:1:3) === 0:3
	@test unify_grades(4, reverse(0:2:4)) === 0:2:4
	@test unify_grades(2, -10:10) === 0:2
	@test unify_grades(2, (2, 0)) === (0, 2)

	@test unify_grades(4, 4, 4) === 4
	@test unify_grades(4, 4, 0, (4, 0)) === (0, 4) # scalar-pseudoscalar subalgebra
	@test unify_grades(4, 4, 2) === 0:2:4 # even subalgebra
	@test unify_grades(4, 4, 1) === 0:4 # full algebra

	@test unify_grades(3, 0:3, 5) === 0:3
	@test unify_grades(3, 2, 0:2:3) === 0:2:3
	@test unify_grades(3, (3, 0), 0) === (0, 3)
	@test unify_grades(3, 0:2:3, 1) === 0:3
end

@testset "grade inference" begin
	for n in 0:4,
		p in [0:n..., 0:n, 0:2:n, n > 0 ? (0, n) : 0],
		q in [0:n..., 0:n, 0:2:n, n > 0 ? (0, n) : 0]

		pq = unify_grades(n, p, q)
		@test p ⊆ pq ⊇ q
		@test length(pq) != 1 || pq isa Integer

		a = resulting_multivector_type(+, Multivector{n,p}, Multivector{n,q})
		@test p ⊆ grade(a) ⊇ q
	end
end

@testset "getindex" begin
	m = Multivector{3,0:3}(1:8)

	@test m[0] === Multivector{3,0}(1:1)
	@test m[1] === Multivector{3,1}(2:4)
	@test m[2:2] === Multivector{3,2}(5:7)
	@test m[(3,)] === Multivector{3,3}(8:8)

	@test m[0:3] == m
	@test m[(0, 3)] == Multivector{3,(0,3)}([1, 8])

	m = zero(Multivector{4,0:2:4,Vector{Int}})

	m[2].comps .= 1:6
	@testset m.comps == [0; 1:6; 0]
end