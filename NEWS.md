# News

## v0.2.0

* Remove `KVector{Sig,K}` and `Multivector{Sig}` in favour of `Multivector{Sig,K}`, where the grade `K` may now be a collection (e.g., range, tuple) in order to represent inhomogeneous multivectors.
This simplifies and generalises the types: only one parametric type is needed, and arbitrary grade combinations (e.g., `0:2:dim` for even multivectors) have efficient representations.

## v0.1.2

* Use `Symbolics.jl` to optimise geometric and derived products by performing computations algebraically and reusing the results.

* Naive implementations of `log`, `sqrt` and trig functions using matrix methods.

* Add duality operations `flipdual`, `hodgedual` and `poincaredual`. Have not defined `dual` because there doesn’t seem to be a single universal notion of multivector duality.

* Rename `Blade`, `Multivector` and `MixedMultivector` to `BasisBlade`, `KVector` and `Multivector`, respectively, to more closely reflect standard terminology.

	- `BasisBlade` is appropriate because only basis blades are representable, not all blades.
	- `KVector` aligns with “k-vector”, common for homogeneous multivectors.
	- `Multivector` is a general term which includes inhomogeneous multivectors.

## v0.1.1

* Add inner and scalar products and left/right contraction products.

* Add `cayleytable` for printing multiplication tables, inspired by [ATell-SoundTheory/CliffordAlgebras.jl](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl).

* Add `vector_repr` and `matrix_repr` for converting between geometric algebra and matrix algebra. (Only the canonical linear representation is supported.)

* Add `Cl(p, q, r)` metric signature, and allow integers for Euclidean spaces.

## v0.1.0
