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
function Base.show(io::IO, T::Type{<:AbstractMultivector{sig,C}}) where {sig,C}
	# important to expose C so T::DataType not ::UnionAll
	name = nameof(T)
	pretty_sig = show_signature(sig)
	params = T.parameters[2:end] # assumes T::DataType
	print(io, "$name{$pretty_sig, $(join(params, ", "))}")
end

show_header(io::IO, a::MixedMultivector) = println(io, typeof(a))
show_header(io::IO, a::HomogeneousMultivector{k}) where k = println(io, "$k-$(typeof(a))")



#= BLADES, MULTIVECTORS, MIXEDMULTIVECTORS =#

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

function show_blade(io::IO, sig, coeff, ublade; styled=false)
	# show coeff with parentheses if necessary
	Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
	print(io, COEFF_BASIS_SEPARATOR)
	isfirst = true
	for bv ∈ each_ublade_bv(ublade)
		label = bvlabel(sig, bv)
		if isfirst
			isfirst = false
		else
			# print(io, length(label) > 1 ? "∧" : BASIS_SEPARATOR)
			print(io, BASIS_SEPARATOR)
		end
		if styled
			printstyled(io, label; bold=true)
		else
			print(io, label)
		end
	end
end

function show_multivector(io::IO, m::Multivector; inline, indent=0, onlynonzero=:auto)
	if onlynonzero == :auto
		onlynonzero = length(m.comps) > 10
	end
	isfirst = true
	for (ublade, coeff) ∈ comps(m)
		if onlynonzero && iszero(coeff) continue end
		if isfirst
			isfirst = false
		else
			print(io, inline ? " + " : "\n")
		end
		print(io, " "^indent)
		show_blade(io, signature(m), coeff, ublade)
	end
	isfirst && print(io, " "^indent, zero(eltype(m)))
end


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
		show_multivector(io, mk; inline=true, onlynonzero=true)
		showparens && print(io, ")")
	end
end


Base.show(io::IO, b::Blade) = show_blade(io, signature(b), b.coeff, b.ublade)
function Base.show(io::IO, ::MIME"text/plain", b::Blade)
	show_header(io, b)
	print(io, " ")
	show(io, b)
end

Base.show(io::IO, m::Multivector) = show_multivector(io, m; inline=true)
function Base.show(io::IO, ::MIME"text/plain", m::Multivector)
	show_header(io, m)
	show_multivector(io, m; inline=false, indent=1)
end

Base.show(io::IO, m::MixedMultivector) = show_mixedmultivector(io, m; inline=true)
function Base.show(io::IO, ::MIME"text/plain", m::MixedMultivector)
	show_header(io, m)
	show_mixedmultivector(io, m; inline=false, indent=1)
end