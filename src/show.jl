DEFAULT_BASIS_SYMBOL = "v"
BASIS_SEPARATOR = ""
COEFF_BASIS_SEPARATOR = " "

const subscript_nums = '₀':'₉'
subscriptnum(n::Integer) = join(subscript_nums[begin + i] for i ∈ reverse(digits(n)))

#= SIGNATURES

Signatures are displayed in shorthand (e.g., (-1,1,1,1) => "⟨-+++⟩") when appearing
in type signatures, but not when displayed by themselves (e.g., with signature(x))
because we don't want it to disguise the actual content.

Subtypes of `AbstractMetricSignature` have a short form defined by `show_signature`
but display in full form when viewed by themselves.

E.g.:
julia> sig = EuclideanSignature(5)
⟨5+⟩ = EuclideanSignature(5)

julia> AbstractMultivector{sig}
AbstractMultivector{⟨5+⟩, C} where C

=#

Base.show(io::IO, sig::AbstractMetricSignature) = print(io, show_signature(sig))
function Base.show(io::IO, ::MIME"text/plain", sig::AbstractMetricSignature)
	print(io, "$sig = ")
	Base.show_default(io, sig)
end

# WARNING: this disguises tuple signatures (1, 1, 1) in type info for the sake of prettiness
function Base.show(io::IO, T::Type{<:AbstractMultivector})
	if isconcretetype(T)
		name = nameof(T)
		pretty_sig = show_signature(signature(T))
		params = T.parameters[2:end]
		print(io, "$name{$pretty_sig, $(join(params, ", "))}")
	else
		invoke(show, Tuple{IO,Type}, io, T) # call original show method
	end
end




#= BLADES, MULTIVECTORS, MIXEDMULTIVECTORS =#

show_header(io::IO, a::MixedMultivector) = println(io, "$(typeof(a)):")
show_header(io::IO, a::HomogeneousMultivector{k}) where k = println(io, "$k-grade $(typeof(a)):")

bvlabel(sig, i) = "$DEFAULT_BASIS_SYMBOL$(subscriptnum(i))"
bvlabel(sig::NamedTuple, i::Integer) = string(keys(sig)[i])
bvlabel(sig, i::Symbol) = string(i)
bvlabel(::MetricSignature{sig}, i) where sig = bvlabel(sig, i)
bvlabel(::OffsetSignature{sig}, i) where sig = bvlabel(sig, i)


# get indices of basis vectors belonging to given ublade
function each_ublade_bv(u::Unsigned)
	bvs = []
	bv = 1
	while u > 0
		if !iszero(u & 1)
			push!(bvs, bv)
		end
		u >>= 1
		bv += 1
	end
	bvs
end
each_ublade_bv(u::Vector) = u


"""
Display unit blade without coefficient.
"""
function show_ublade(io::IO, sig, ublade)
	isfirst = true
	for bv ∈ each_ublade_bv(ublade)
		label = bvlabel(sig, bv)
		if isfirst
			isfirst = false
		else
		print(io, BASIS_SEPARATOR)
		end
		if get(io, :color, false)
			printstyled(io, label; bold=true)
		else
			print(io, label)
		end
	end
end


"""
Display blade with parentheses surrounding coefficient if necessary.
"""
function show_blade(io::IO, sig, coeff, ublade)
	Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
	ublade_grade(ublade) == 0 && return
	print(io, COEFF_BASIS_SEPARATOR)
	show_ublade(io, sig, ublade)
end

sorted_blades(a::Blade) = blades(a)
sorted_blades(a::CompositeMultivector{sig,<:AbstractVector}) where sig = blades(a)
sorted_blades(a::CompositeMultivector{sig,<:AbstractDict}) where sig = sort(blades(a))

function show_multivector_inline(io::IO, m::Multivector; onlynonzero=:auto)
	isfirst = true
	for u ∈ sorted_blades(m)
		onlynonzero == true && iszero(u) && continue
		isfirst ? isfirst = false : print(io, " + ")
		show(io, u)
	end
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
function show_multivector(io::IO, m::Multivector; indent=0, onlynonzero=:auto)
	if iszero(m)
		print(io, " "^indent, zero(eltype(m)))
		return
	end
	mcomps = [(u.ublade, u.coeff) for u ∈ blades(m)]
	alignments = [Base.alignment(io, coeff) for (_, coeff) ∈ mcomps]
	L = maximum(first.(alignments))
	R = maximum(last.(alignments))
	isfirst = true
	for ((u, coeff), (l, r)) ∈ zip(mcomps, alignments)
		isfirst ? isfirst = false : println(io)
		print(io, " "^(L - l + indent))
		Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
		print(io, " "^(R - r), " ")
		show_ublade(io, signature(m), u)
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
	for k ∈ 0:maximum(grade(m))
		mk = grade(m, k)
		if iszero(mk) continue end
		if firstline
			firstline = false
		else
			print(io, inline ? " + " : "\n")
		end
		print(io, " "^indent)
		showparens = inline && (count(!iszero, mk.comps) > 1)
		showparens && print(io, "(")
		show_multivector_inline(io, mk; onlynonzero=true)
		showparens && print(io, ")")
	end
end


Base.show(io::IO, b::Blade) = show_blade(io, signature(b), b.coeff, b.ublade)
function Base.show(io::IO, ::MIME"text/plain", b::Blade)
	show_header(io, b)
	print(io, " ")
	show(io, b)
end

Base.show(io::IO, m::Multivector) = show_multivector_inline(io, m)
function Base.show(io::IO, ::MIME"text/plain", m::Multivector)
	show_header(io, m)
	show_multivector(io, m; indent=1)
end

Base.show(io::IO, m::MixedMultivector) = show_mixedmultivector(io, m; inline=true)
function Base.show(io::IO, ::MIME"text/plain", m::MixedMultivector)
	show_header(io, m)
	show_mixedmultivector(io, m; inline=false, indent=1)
end