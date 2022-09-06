@testset "show methods" begin
	io = IOBuffer()
	for a in [
		Blade{(1,1)}(0b00 => 42),
		Blade{(1,)}(0b1 => 1.1),
		Blade{(1,1,1)}(0b111 => true),
	], b in [a, Multivector(a), MixedMultivector(a)]
		@test isnothing(show(io, b))
		@test isnothing(show(io, MIME("text/plain"), b))
		@test repr(b) isa String
	end
end