"""
Display unit basis blade without coefficient.
"""
function show_basis_blade(io::IO, sig, bits)
	# blades with higher grade than dimension are always zero,
	#  and do not need to have any basis printed
	grade(bits) > dimension(sig) && return
	indices = bits_to_indices(bits)
	printstyled(io, basis_blade_label(sig, indices); bold=true)
end
show_basis_blade(io::IO, b::Blade) = show_basis_blade(io, signature(b), bitsof(b))

plaintext_repr(io, a) = repr(a; context=IOContext(io, :color=>false))

function Base.alignment(io::IO, b::Blade)
	(l, r) = Base.alignment(io, b.coeff)
	(l, length(plaintext_repr(io, b)) - l)
end

function Base.alignment(io::IO, b::CompositeMultivector)
	(length(plaintext_repr(io, b)), 1)
end

COEFF_BASIS_SEPARATOR = " "

"""
Display blade with parentheses surrounding coefficient if necessary.
"""
function show_blade(io::IO, b::Blade; compact=false)
	Base.show_unquoted(io, b.coeff, 0, Base.operator_precedence(:*))
	grade(b) == 0 && return
	compact || print(io, COEFF_BASIS_SEPARATOR)
	(iszero(b) && compact) && print(io, '*')
	show_basis_blade(io, b)
end


"""
Display a multivector as a column of blades, with coefficients aligned using
the native alignment mechanism, and blades basis aligned.

```
julia> 1e3x + y + 1e-3z
1-Multivector{⟨x+,y+,z+⟩, Vector{Float64}, 1}
 1000.0   x
    1.0   y
    0.001 z
```
"""
function show_multivector(io::IO, m::Multivector; indent=0)
	if iszero(m)
		print(io, " "^indent, zero(eltype(m)))
		return
	end
	mcomps = [(bitsof(u), u.coeff) for u ∈ blades(m)]
	alignments = [Base.alignment(io, coeff) for (_, coeff) ∈ mcomps]
	L = maximum(first.(alignments))
	R = maximum(last.(alignments))
	isfirst = true
	for ((u, coeff), (l, r)) ∈ zip(mcomps, alignments)
		isfirst ? isfirst = false : println(io)
		print(io, " "^(L - l + indent))
		Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
		print(io, " "^(R - r), " ")
		show_basis_blade(io, signature(m), u)
	end
end


sorted_blades(a::Blade) = blades(a)
sorted_blades(a::CompositeMultivector{<:AbstractVector}) = blades(a)
sorted_blades(a::CompositeMultivector{<:AbstractDict}) = sort(blades(a), by=(b -> bits(b)))

function show_multivector_inline(io::IO, m::Multivector; compact=false)
	if iszero(m)
		print(io, zero(eltype(m)))
		return
	end
	isfirst = true
	for u ∈ sorted_blades(m)
		compact && iszero(u) && continue
		isfirst ? isfirst = false : print(io, " + ")
		show_blade(io, u; compact)
	end
end





"""
Display inhomogeneous `MixedMultivector` with each grade on a new line.
"""
function show_mixedmultivector(io::IO, m::MixedMultivector; inline, indent=0)
	firstline = true
	if iszero(m)
		print(io, " "^indent, zero(eltype(m)))
		return
	end
	for k ∈ 0:dimension(m)
		mk = grade(m, k)
		if iszero(mk) continue end
		if firstline
			firstline = false
		else
			print(io, inline ? " + " : "\n")
		end
		print(io, " "^indent)
		showparens = inline && (length(blades(mk)) > 1)
		showparens && print(io, "(")
		show_multivector_inline(io, mk; compact=inline)
		showparens && print(io, ")")
	end
end



show_header(io::IO, a::MixedMultivector) = println(io, "$(typeof(a)):")
show_header(io::IO, a::HomogeneousMultivector{sig,k}) where {sig,k} = println(io, "Grade-$k $(typeof(a)):")


Base.show(io::IO, b::Blade) = show_blade(io, b; compact=true)
function Base.show(io::IO, ::MIME"text/plain", b::Blade)
	show_header(io, b)
	print(io, " ")
	show_blade(io, b)
end

Base.show(io::IO, m::Multivector) = show_multivector_inline(io, m; compact=true)
function Base.show(io::IO, ::MIME"text/plain", m::Multivector)
	show_header(io, m)
	show_multivector(io, m; indent=1)
end

Base.show(io::IO, m::MixedMultivector) = show_mixedmultivector(io, m; inline=true)
function Base.show(io::IO, ::MIME"text/plain", m::MixedMultivector)
	show_header(io, m)
	show_mixedmultivector(io, m; inline=false, indent=1)
end