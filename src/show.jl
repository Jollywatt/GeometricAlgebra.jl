
#= Blades =#

# make things line up right when printing in arrays
function Base.alignment(io::IO, @nospecialize(b::BasisBlade))
	(l, r) = Base.alignment(io, b.coeff)
	(l, length(sprint(show, b)) - l)
end

"""
Display blade with parentheses surrounding coefficient if necessary.
	
```jldoctest
julia> GeometricAlgebra.show_blade(stdout, BasisBlade{(x=1,)}(0b1 => 1 + im))
(1+1im) x
```
"""
function show_blade(io::IO, @nospecialize(b::BasisBlade); compact=false)
	subio = IOContext(io, :compact => true)
	if compact && isnumberzero(b.coeff)
		print(io, "0")
		return
	elseif compact && isnumberone(b.coeff)
		isscalar(b) && print(io, "1")
	elseif compact && isnumberone(-b.coeff) 
		print(io, isscalar(b) ? "-1" : "-")
	else
		Base.show_unquoted(subio, b.coeff, 0, Base.operator_precedence(:*))
	end
	grade(b) == 0 && return
	# blades with higher grade than dimension are always zero,
	#  and should not have basis blade printed
	0 < grade(b) <= dimension(b) || return
	compact || print(io, " ") # coefficient–basis separator
	show_basis_blade(io, signature(b), bits_to_indices(bitsof(b)))
end



#= CompositeMultivectors =#

"""
Display homogeneous multivector components as a column of blades,
with coefficients and blades aligned using the native alignment mechanism.

```jldoctest
julia> a = KVector{(1,1,1),1}([1e3, 1, 1e-3]);

julia> GeometricAlgebra.show_multivector(stdout, a)
1000.0   v1
   1.0   v2
   0.001 v3
```
"""
function show_multivector(io::IO, @nospecialize(a); indent=0, showzeros=true)
	# TODO: showzeros argument
	iszero(a) && return print(io, " "^indent, numberzero(eltype(a)))

	comps = zip(bitsof(a), a.comps)
	if !showzeros
		comps = Iterators.filter(!iszero∘last, comps)
	end

	comps = collect(comps)

	alignments = Base.alignment.(Ref(io), last.(comps))
	L = maximum(first.(alignments))
	R = maximum(last.(alignments))

	firstline = true
	for ((bits, coeff), (l, r)) ∈ zip(comps, alignments)
		firstline || println(io)
		print(io, " "^(L - l + indent))
		Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
		print(io, " "^(R - r), " ")
		show_basis_blade(io, signature(a), bits_to_indices(bits))
		firstline = false
	end
end


function show_multivector_inline(io::IO, @nospecialize(a); compact=false, showzeros=false)
	if (!showzeros || compact) && iszero(a)
		print(io, numberzero(eltype(a)))
		return
	end
	isfirst = true
	for (bits, coeff) in zip(bitsof(a), a.comps)
		!showzeros && isnumberzero(coeff) && continue
		isfirst ? isfirst = false : print(io, " + ")
		show_blade(io, BasisBlade{signature(a)}(bits => coeff); compact)
	end
end


"""
Display an inhomogeneous `Multivector` with each grade on a new line.
"""
function show_mixedmultivector(io::IO, @nospecialize(a); inline, indent=0, showzeros=false)
	if iszero(a)
		print(io, " "^indent, numberzero(eltype(a)))
		return
	end
	firstline = true
	for k ∈ grade(a)
		ak = grade(a, k)

		!showzeros && iszero(ak) && continue
		firstline || print(io, inline ? " + " : "\n")

		if firstline || !inline
			print(io, " "^indent)
		end

		showparens = inline && (0 < k < dimension(a))
		showparens && print(io, "(")
		show_multivector_inline(io, ak; compact=inline, showzeros)
		showparens && print(io, ")")

		firstline = false
	end
end



function show_header(io::IO, @nospecialize(a::BasisBlade))
	show(io, typeof(a))
	println(io, ":")
end
function show_header(io::IO, @nospecialize(a::Multivector))
	print(io, ncomponents(a), "-component ")
	show(io, typeof(a))
	println(io, ":")
end

Base.show(io::IO, @nospecialize(a::BasisBlade)) = show_blade(io, a; compact=true)
function Base.show(io::IO, ::MIME"text/plain", @nospecialize(a::BasisBlade))
	show_header(io, a)
	print(io, " ")
	show_blade(io, a)
end

Base.show(io::IO, @nospecialize(a::Multivector)) = show_mixedmultivector(io, a; inline=true)
function Base.show(io::IO, ::MIME"text/plain", @nospecialize(a::Multivector))
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
