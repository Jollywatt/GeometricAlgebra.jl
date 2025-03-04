using GeometricAlgebra.MiniCAS

macro ident(expr)
	@assert expr isa Expr
	@assert expr.head == :call
	op, a, b = expr.args
	@assert op == :(==)
	quote
		@test $a == $b
		@test isequal($a, $b)
		@test hash($a) == hash($b)
	end |> esc
end

@testset "products" begin
	x, y, z = ProductNode.([:x, :y, :z] .=> 1)

	@ident x*y == y*x
	@ident x*x == x^2
	@ident inv(x) == x/x^2
	@ident (x*y)/(x^3*z) == y/z*inv(x)^2
	@ident x/x == one(x)
end


@testset "sums" begin
	x, y, z = ProductNode.([:x, :y, :z] .=> 1)

	@ident 3x == x*2 + x
	@ident x/2 == 0.5x
	@ident 2\x/3 == inv(6)x

	@ident x + y == y + x
	@ident x + x == 2x
	@ident -x == x - 2x

	@ident x - x == y - y

	@ident x*y == y*x
	@ident (x + y)*z == x*z + y*z

	@ident 0y*(x + z) == 0z

	@ident (x + y)*(x - y) == x^2 - y^2
	@ident (x + y)^6/(y + x)^8 == inv(x + z + y - z)^2

	@ident x - 1 == y\(x*y - y)

	a = x + 2y + 3z
	@ident a/x == 1 + 2y/x + x\3z

end

@testset "factor" begin
	x, y, z = ProductNode.([:x, :y, :z] .=> 1)

	a = x*y - 2
	s, a = MiniCAS.fixmul!!(a)
	@test factor(x^2*y - 2x) == SumNode(ProductNode(a => 1, :x => 1) => s)

	@ident factor(4a/-8a*x) == -.5x
	@ident factor(2a*x*y/a) == 2y*x

	@test isone(factor(x*y/a - 2/a))
	@test isone(factor(20/a - 10x*y/a)/-10)

	@ident inv(-2a) == -inv(a)/2

end