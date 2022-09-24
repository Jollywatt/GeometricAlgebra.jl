
#= Blades =#

# make things line up right when printing in arrays
function Base.alignment(io::IO, b::Blade)
	(l, r) = Base.alignment(io, b.coeff)
	(l, length(sprint(show, b)) - l)
end

"""
Display blade with parentheses surrounding coefficient if necessary.
	
```jldoctest
julia> GeometricAlgebra.show_blade(stdout, Blade{(x=1,)}(0b1 => 1 + im))
(1+1im) x
```
"""
function show_blade(io::IO, b::Blade; compact=false)
	subio = IOContext(io, :compact => true)
	if compact && isrealzero(b.coeff) print(io, "0")
	elseif compact && isrealone(b.coeff) isscalar(b) && print(io, "1")
	elseif compact && isrealone(-b.coeff) print(io, isscalar(b) ? "-1" : "-")
	else
		Base.show_unquoted(subio, b.coeff, 0, Base.operator_precedence(:*))
	end
	grade(b) == 0 && return
	compact || print(io, " ") # coefficient–basis separator
	# blades with higher grade than dimension are always zero,
	#  and should not have basis blade printed
	0 < grade(b) <= dimension(b) || return
	show_basis_blade(io, signature(b), bits_to_indices(bitsof(b)))
end



#= CompositeMultivectors =#

"""
Display homogeneous multivector components as a column of blades,
with coefficients and blades aligned using the native alignment mechanism.

```jldoctest
julia> a = Multivector{(1,1,1),1}([1e3, 1, 1e-3]);

julia> GeometricAlgebra.show_multivector(stdout, a)
1000.0   v1
   1.0   v2
   0.001 v3
```
"""
function show_multivector(io::IO, a::Multivector; indent=0)
	iszero(a) && return print(io, " "^indent, realzero(eltype(a)))

	alignments = Base.alignment.(Ref(io), a.components)
	L = maximum(first.(alignments))
	R = maximum(last.(alignments))

	for (i, (bits, (l, r))) ∈ enumerate(zip(bits_of_grade(grade(a), dimension(a)), alignments))
		i > 1 && println(io)
		print(io, " "^(L - l + indent))
		Base.show_unquoted(io, a.components[i], 0, Base.operator_precedence(:*))
		print(io, " "^(R - r), " ")
		show_basis_blade(io, signature(a), bits_to_indices(bits))
	end
end


function show_multivector_inline(io::IO, a::Multivector; compact=false, showzeros=false)
	if iszero(a)
		print(io, realzero(eltype(a)))
		return
	end
	isfirst = true
	for (bits, coeff) in zip(bits_of_grade(grade(a), dimension(a)), a.components)
		(!showzeros || compact) && isrealzero(coeff) && continue
		isfirst ? isfirst = false : print(io, " + ")
		show_blade(io, Blade{signature(a)}(bits => coeff); compact)
	end
end


"""
Display an inhomogeneous `MixedMultivector` with each grade on a new line.
"""
function show_mixedmultivector(io::IO, a::MixedMultivector; inline, indent=0)
	firstline = true
	if iszero(a)
		print(io, " "^indent, realzero(eltype(a)))
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



function show_header(io::IO, a::Blade)
	show(io, typeof(a))
	println(io, ":")
end
function show_header(io::IO, a::CompositeMultivector)
	print(io, length(a.components), "-component ")
	show(io, typeof(a))
	println(io, ":")
end

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



#= Pretty-printing =#


# display object with given show function along with unabbreviated form if different
function show_pretty(io, showfn, x)
	compact = sprint(showfn, x, context = :compact => true)
	full = sprint(show, x, context = :compact => false)
	showfn(io, x)
	compact == full && return
	printstyled(io, " (pretty-printed ", full, ")", color = :light_black)
end



#= Metric Signatures =#


struct SignatureDisplayWrapper{Sig} end
Base.show(io::IO, ::SignatureDisplayWrapper{Sig}) where {Sig} = show_signature(io, Sig)

# use `show_signature` to pretty-print the first type parameter of multivector types
# by wrapping `Sig` in `SignatureDisplayWrapper{Sig}()`
function pretty_print_signature(io, T::Type)
	if get(io, :compact, true)::Bool && Base.unwrap_unionall(T) isa DataType

		# unwrap UnionAll and remember TypeVars
		typevars = TypeVar[]
		T = T
		while T isa UnionAll
			push!(typevars, T.var)
			T = T.body
		end

		sig, other_parameters... = T.parameters
		# wrap sig type parameter in pretty-printer
		sig = sig isa TypeVar ? sig : SignatureDisplayWrapper{sig}()
		T = T.name.wrapper{sig, other_parameters...}

		# rewrap UnionAll in TypeVars
		for p in reverse(typevars)
			T = UnionAll(p, T)
		end
	end

	invoke(show, Tuple{IO,Type}, io, T) # show with default method
end




# when type is displayed on its own, show in full to avoid confusing the user
Base.show(io::IO, T::Type{<:AbstractMultivector}) = pretty_print_signature(io, T)
function Base.show(io::IO, ::MIME"text/plain", T::Type{<:AbstractMultivector})
	show_pretty(io, show, T)
end
