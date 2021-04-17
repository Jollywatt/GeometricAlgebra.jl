basis_vectors(sig) = [Blade{sig}(1, i) for i ∈ Iterators.take(FixedGradeBits(1), dimension(sig))]
basis(sig) = basis_vectors(sig)
basis(osig::OffsetSignature{sig,indices}) where {sig,indices} = OffsetArray(basis_vectors(osig), indices)
basis(dim::Integer) = basis(Tuple(fill(1, dim)))