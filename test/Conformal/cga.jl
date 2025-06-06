using GeometricAlgebra.Conformal

@testset "blade standardisation" begin

	@testset for dim in 0:4
		n0, noo = nullbasis(dim)

		@testset for k in 0:dim
			for _ in 1:10
				E = k > 0 ? wedge(randn(Multivector{dim,1}, k)...) : randn()
				p = randn(Multivector{dim,1})
				r² = randn()

				dir = E∧noo
				flat = translate(p, n0∧E∧noo)
				dualflat = translate(p, E)
				round = translate(p, (n0 + 2\r²*noo)∧E)

				dirblade = standardform(dir)
				flatblade = standardform(flat)
				dualflatblade = standardform(dualflat)
				roundblade = standardform(round)

				@test dirblade isa Conformal.DirectionBlade
				@test flatblade isa Conformal.FlatBlade
				@test dualflatblade isa Conformal.DualFlatBlade
				@test roundblade isa Conformal.RoundBlade

				@test Multivector(dirblade) ≈ dir
				@test Multivector(flatblade) ≈ flat
				@test Multivector(dualflatblade) ≈ dualflat
				@test Multivector(roundblade) ≈ round
			end
		end
	end
end
