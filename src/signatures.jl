dimension(sig) = length(sig)

basis_labels(sig::Tuple, indices=1:dimension(sig)) = ["v$i" for i ∈ indices]
basis_labels(sig::NamedTuple{labels}) where labels = labels

basis_blade_label(sig::Tuple, indices) = "v"*join(string.(indices))
basis_blade_label(sig, indices) = join(basis_labels(sig)[indices])


struct OffsetSignature{sig,indices} end
OffsetSignature(sig,indices) = OffsetSignature{sig,indices}()
Base.getindex(::OffsetSignature{sig}, args...) where sig = getindex(sig, args...)
Base.length(::OffsetSignature{sig}) where sig = length(sig)
basis_labels(::OffsetSignature{sig,indices}) where {sig,indices} = basis_labels(sig, indices)


"""
Convert human-readable basis vector index or symbol into 1-based integer index
"""
normalize_bv_index(sig, i) = i
normalize_bv_index(sig, i::Symbol) = findfirst(==(i), basis_labels(sig))
function normalize_bv_index(sig::Union{Tuple,NamedTuple}, i::Integer)
	r = 1:dimension(sig)
	i ∈ r ? i : error("index $i is outside range $r for signature $sig")
end
function normalize_bv_index(osig::OffsetSignature{sig,indices}, ioffset::Integer) where {sig,indices}
	i = findfirst(==(ioffset), indices)
	isnothing(i) && error("index $ioffset is outside range $indices for $osig")
	i
end


"""
Pretty-print metric signature in short, non-parseable form such as ⟨-+++⟩
"""
show_signature(sig) = repr(sig) # fallback
show_signature(sig::Tuple) = "⟨$(join(map(s -> get(Dict(+1=>"+", -1=>"-"), s, s), sig)))⟩"
function show_signature(sig::NamedTuple)
	items = map(zip(keys(sig), sig)) do (label, square)
		s = get(Dict(+1=>"+", -1=>"-"), square, square)
		"$label$s"
	end
	"⟨$(join(items, ","))⟩"
end
function show_signature(::OffsetSignature{sig,indices}) where {sig,indices}
	"$(show_signature(sig)) with indices $indices"
end