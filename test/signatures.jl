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