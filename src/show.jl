#= PRETTY-PRINTED TYPES

Magic is used to pretty-print the metric signature as it appears in a type:
`MixedMultivector{(1, 1, 1), S}` displays as `MixedMultivector{⟨+++⟩, S} where S`
=#



"""
	ShowWrapper

Wrapper around an all-bits object which displays as if the object was printed directly.

Examples
===
```jldoctest
julia> GeometricAlgebra.ShowWrapper("⟨+++⟩")
⟨+++⟩

julia> typeof(ans)
GeometricAlgebra.ShowWrapper{Symbol("⟨+++⟩")}
```
"""
struct ShowWrapper{object} end
ShowWrapper(a) = ShowWrapper{Symbol(a)}()
ShowWrapper(a::Symbol) = ShowWrapper(repr(a))

Base.show(io::IO, sw::ShowWrapper{object}) where object = print(io, object)


"""
	pretty_print_type_parameters(io, T, fns...)

Print a `Type` (`DataType` or `UnionAll`) with each parameter pretty-printed
if it is specified. The `i`th type parameter is displayed as what `fns[i]`
returns when called with the parameter.

Examples
===
```jldoctest; setup = :( using GeometricAlgebra: pretty_print_type_parameters )
julia> pretty_print_type_parameters(stdout, Array{T,2} where T,
       	identity,
       	n -> "\$n-dimensional")
Array{T, 2-dimensional} where T
```
"""
function pretty_print_type_parameters(io::IO, T::Type, fns...)
	typevars = TypeVar[]
	while T isa UnionAll # unwrap UnionAll and remember TypeVars
		push!(typevars, T.var)
		T = T.body
	end
	parameters = [p isa TypeVar ? p : ShowWrapper(fn(p)) for (p, fn) ∈ zip(T.parameters, fns)]
	S = T.name.wrapper{parameters...}
	for tv in reverse(typevars) # rewrap into UnionAll
		S = UnionAll(tv, S)
	end
	invoke(show, Tuple{IO,Type}, io, S) # show with default method
end



# pretty-print signature and display bits parameter as short binary string (instead of long hexadecimal number)
function Base.show(io::IO, T::Type{<:Blade})
	dim = isconcretetype(T) ? dimension(T) : 0
	pretty_print_type_parameters(io, T,
		show_signature,
		identity,
		function(bits)
			bits isa Unsigned || return bits
			"0b"*string(bits, pad=dim, base=2)
		end,
		identity
	)
end
Base.show(io::IO, T::Type{<:Multivector}) = pretty_print_type_parameters(io, T, show_signature, identity, identity)
Base.show(io::IO, T::Type{<:MixedMultivector}) = pretty_print_type_parameters(io, T, show_signature, identity)

# when type is displayed on its own, be explicit about pretty-printing to avoid confusing the user
function Base.show(io::IO, ::MIME"text/plain", T::Type{<:AbstractMultivector})
	show(io, T)
	Base.print_without_params(T) && return
	println(io)
	printstyled(io, "(pretty-printed ")
	invoke(show, Tuple{IO,Type}, io, T) # show with default method
	printstyled(io, ")")
end

# when signatures are displayed on their own, show both the shorthand and the full expression
function Base.show(io::IO, ::MIME"text/plain", sig::MetricSignature)
	println(io, show_signature(sig))
	printstyled(io, "(pretty-printed ")
	show(io, sig) # show with default method
	printstyled(io, ")")
end



#= BLADE DISPLAY METHODS =#

# show basis blade only without scalar coefficient
function show_basis_blade(io::IO, sig, bits)
	# blades with higher grade than dimension are always zero,
	#  and do not need to have any basis printed
	0 < grade(bits) <= dimension(sig) || return
	label = basis_blade_label(sig, bits_to_indices(bits))
	printstyled(io, label; bold=true)
end
show_basis_blade(io::IO, b::Blade) = show_basis_blade(io, signature(b), bitsof(b))

plaintext_repr(io, a) = repr(a; context=IOContext(io, :color=>false))

# make things line up right when printing in arrays
function Base.alignment(io::IO, b::Blade)
	(l, r) = Base.alignment(io, b.coeff)
	(l, length(plaintext_repr(io, b)) - l)
end
function Base.alignment(io::IO, b::CompositeMultivector)
	(0, length(plaintext_repr(io, b)))
end

SHOW_IOCONTEXT = Dict(:compact => true)

"""
Display blade with parentheses surrounding coefficient if necessary.
	
```jldoctest
julia> GeometricAlgebra.show_blade(stdout, Blade{(x=1,)}(1 + im, 0b1))
(1+1im) x
```
"""
function show_blade(io::IO, b::Blade; compact=false)
	subio = IOContext(io, SHOW_IOCONTEXT...)
	Base.show_unquoted(subio, b.coeff, 0, Base.operator_precedence(:*))
	grade(b) == 0 && return
	compact || print(io, " ") # coefficient--basis separator
	(iszero(b) && compact) && print(io, '*') # write `0*x` instead of `0x`, etc
	show_basis_blade(io, b)
end



#= COMPOSITE MULTIVECTOR DISPLAY METHODS =#

"""
Display a multivector as a column of blades, with coefficients aligned using
the native alignment mechanism, and blades basis aligned.

```jldoctest; setup = :( (x, y, z) = basis((x=1, y=1, z=1)) )
julia> GeometricAlgebra.show_multivector(stdout, 1e3x + y + 1e-3z)
1000.0   x
   1.0   y
   0.001 z
```
"""
function show_multivector(io::IO, a::Multivector; indent=0)
	iszero(a) && return print(io, " "^indent, zero(eltype(a)))

	comps = [(bitsof(u), u.coeff) for u ∈ blades(a)]
	alignments = [Base.alignment(io, coeff) for (_, coeff) ∈ comps]
	L = maximum(first.(alignments))
	R = maximum(last.(alignments))
	isfirst = true
	for ((bits, coeff), (l, r)) ∈ zip(comps, alignments)
		isfirst ? isfirst = false : println(io)
		print(io, " "^(L - l + indent))
		Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
		print(io, " "^(R - r), " ")
		show_basis_blade(io, signature(a), bits)
	end
end


function show_multivector_inline(io::IO, a::Multivector; compact=false, showzeros=false)
	if iszero(a)
		print(io, zero(eltype(a)))
		return
	end
	isfirst = true
	for b ∈ blades(a)
		(!showzeros || compact) && iszero(b) && continue
		isfirst ? isfirst = false : print(io, " + ")
		show_blade(io, b; compact)
	end
end


"""
Display inhomogeneous `MixedMultivector` with each grade on a new line.
"""
function show_mixedmultivector(io::IO, a::MixedMultivector; inline, indent=0)
	firstline = true
	if iszero(a)
		print(io, " "^indent, zero(eltype(a)))
		return
	end
	for k ∈ 0:dimension(a)
		ak = grade(a, k)
		if iszero(ak) continue end
		if firstline
			firstline = false
		else
			print(io, inline ? " + " : "\n")
		end
		print(io, " "^indent)
		showparens = inline && (0 < k < dimension(a))
		showparens && print(io, "(")
		show_multivector_inline(io, ak; compact=inline)
		showparens && print(io, ")")
	end
end



show_header(io::IO, a::MixedMultivector) = println(io, "$(typeof(a)):")
show_header(io::IO, a::HomogeneousMultivector{sig,k}) where {sig,k} = println(io, "Grade-$k $(typeof(a)):")


Base.show(io::IO, a::Blade) = show_blade(io, a; compact=true)
function Base.show(io::IO, ::MIME"text/plain", a::Blade)
	show_header(io, a)
	print(io, " ")
	show_blade(io, a)
end

Base.show(io::IO, a::Multivector) = show_multivector_inline(io, a; compact=true)
function Base.show(io::IO, ::MIME"text/plain", a::Multivector)
	show_header(io, a)
	show_multivector(io, a; indent=1)
end

Base.show(io::IO, a::MixedMultivector) = show_mixedmultivector(io, a; inline=true)
function Base.show(io::IO, ::MIME"text/plain", a::MixedMultivector)
	show_header(io, a)
	show_mixedmultivector(io, a; inline=false, indent=1)
end