using GeometricAlgebra

v1, v2 = basis((1, 1))
v12 = v1*v2

x, y = basis((x=1, y=1))
xy = x*y

@testset "repr parseability" begin
	isparseable(a) = eval(Meta.parse(repr(a))) == a
	
	@test isparseable(4v1)
	@test isparseable(0v1)
	@test isparseable((1 + im)*v1)
	@test isparseable(v1 - 2v2)
	@test isparseable(4 + v12)

	@test isparseable(10x)
	@test isparseable(0*x)
	@test isparseable(0.0x)
	@test isparseable(1 + 2x + 3xy)
end

# just test against errors
@test isnothing(show(stdout, MIME("text/plain"), x))
@test isnothing(show(stdout, MIME("text/plain"), x + y))
@test isnothing(show(stdout, MIME("text/plain"), 1 + y))