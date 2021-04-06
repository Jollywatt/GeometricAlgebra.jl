dimension(sig) = length(sig)

basis_labels(sig::Tuple) = ["v$i" for i ∈ 1:N]
basis_labels(sig::NamedTuple{labels}) where labels = string.(labels)

basis_blade_label(sig::Tuple, indices) = "v"*join(string.(indices))
basis_blade_label(sig, indices) = join(basis_labels(sig)[indices])
