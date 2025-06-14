mutable struct DisplayOptions
	inline::Bool
	groupgrades::Union{Bool,Nothing}
	showzeros::Union{Bool,Nothing}
	eps::Float64
	parseable::Bool
	compact::Bool
end

const display_options = DisplayOptions(
	false,
	nothing,
	nothing,
	0,
	false,
	false,
)

#= BasisBlade =#

# make things line up right when printing in arrays
function Base.alignment(io::IO, @nospecialize(b::BasisBlade))
	style = get_basis_display_style(signature(b))
	coeff = b.coeff*(-1)^basis_blade_parity(style, b.bits)
	coeff_width = length(sprint(show, coeff; context=io))
	basis_width = length(sprint(show_basis_blade, style, b.bits; context=io))
	(coeff_width, basis_width)
end

"""
Display blade with parentheses surrounding coefficient if necessary.

# Example
```jldoctest
julia> GeometricAlgebra.show_blade(stdout, BasisBlade{(x=1,)}(1 + im, 0b1))
(1 + 1im) x
```
"""
function show_blade(io::IO, @nospecialize(a::BasisBlade);
                    compact=get(io, :compact, false),
                    parseable=false,
                    basis_display_style=get_basis_display_style(signature(a)))
	if parseable
		@static if VERSION ≥ v"1.7"
			type = Base.typeinfo_implicit(typeof(a.coeff)) ? constructor(a) : typeof(a)
		else
			type = typeof(a)
		end
		print(io, type, "(", a.coeff, ", ")
		niceshow(io, dimension(a), a.bits)
		print(io, ")")
	else
		parity = basis_blade_parity(basis_display_style, a.bits)
		coeff = (-1)^parity*a.coeff

		if compact && isnumberzero(coeff)
			print(io, "0")
			return
		elseif compact && isnumberone(coeff)
			isscalar(a) && print(io, "1")
		elseif compact && isnumberone(-coeff) 
			print(io, isscalar(a) ? "-1" : "-")
		else
			Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))
		end
		# blades with higher grade than dimension are always zero,
		#  and should not have basis blade printed
		0 < grade(a) <= dimension(a) || return
		compact || print(io, " ") # coefficient–basis separator

		show_basis_blade(io, basis_display_style, a.bits)
	end
end



#= Multivector =#

issmall(x::Number, eps) = abs(x) <= eps
issmall(x, eps) = false
issmall(x::AbstractArray, eps) = all(Base.Fix2(issmall, eps), x)
issmall(x::Multivector, eps) = issmall(x.comps, eps)


function show_multivector_row(io::IO, @nospecialize(a);
	                           indent=0,
	                           compact=false,
	                           showzeros=false,
	                           eps=0,
	                           basis_display_style=get_basis_display_style(signature(a)))
	io = IOContext(io, :compact=>true)
	if (!showzeros || compact) && iszero(a)
		print(io, numberzero(eltype(a)))
		return
	end
	isfirst = true
	for bits ∈ componentbits(basis_display_style, grade(a))
		coeff = a.comps[componentindex(a, bits)]
		!showzeros && issmall(coeff, eps) && continue
		isfirst ? print(io, " "^indent) : print(io, " + ")
		show_blade(io, BasisBlade{signature(a)}(coeff, bits); compact, basis_display_style)
		isfirst = false
	end
end

function show_multivector_col(io::IO, @nospecialize(a);
	                           indent=0,
	                           showzeros=true,
	                           compact=false,
	                           eps=0,
	                           basis_display_style=get_basis_display_style(signature(a)))
	!showzeros && iszero(a) && return print(io, " "^indent, numberzero(eltype(a)))

	bits = componentbits(basis_display_style, grade(a))
	comps = [
		(-1)^basis_blade_parity(basis_display_style, b)*a.comps[componentindex(a, b)] => b
	for b in bits]

	!showzeros && filter!(!Base.Fix2(issmall, eps)∘first, comps)

	isempty(comps) && return print(io, " "^indent, a.comps)

	alignments = Base.alignment.(Ref(io), first.(comps))
	L = maximum(first.(alignments))
	R = maximum(last.(alignments))

	firstline = true
	for ((coeff, b), (l, r)) ∈ zip(comps, alignments)
		firstline || println(io)
		firstline = false

		print(io, " "^(L - l + indent))
		Base.show_unquoted(io, coeff, 0, Base.operator_precedence(:*))

		iszero(b) && continue
		print(io, " "^(R - r), compact ? "" : " ")
		show_basis_blade(io, basis_display_style, b)
	end
end


"""
	show_multivector(io, a; kwargs...)

Display multivector components in a column or inline, optionally grouping by grade.

# Keyword arguments

- `inline::Bool`: print on one line (default `true`).
- `groupgrades::Bool`: visually group components by grade (default `true`).
   If inline, prints parentheses around parts of each grade; if multiline, prints
   each grade on its own line.
- `showzeros::Bool`: whether to omit zero components from display.
- `eps::Real`: supprxess components with a magnitude smaller than a threshold (default `0`).
- `indent::Integer`: indentation width.
- `parseable::Bool`: use parseable style (used by `repr`) instead of human-readable style.
- `compact::Bool`: omit unnecessary spaces, coefficients of unity, etc (default: `false`).
- `basis_display_style::BasisDisplayStyle`: style to use to print basis blades.

# Examples
```jldoctest
julia> a = Multivector{2,0:2}((1:4) .^ 2);

julia> GeometricAlgebra.show_multivector(stdout, a; inline=false, groupgrades=false)
 1
 4 v1
 9 v2
16 v12

julia> GeometricAlgebra.show_multivector(stdout, a; inline=false, groupgrades=true)
1
4 v1 + 9 v2
16 v12

julia> GeometricAlgebra.show_multivector(stdout, a; inline=true, groupgrades=true)
(1) + (4 v1 + 9 v2) + (16 v12)

julia> GeometricAlgebra.show_multivector(stdout, a; inline=true, groupgrades=false)
1 + 4 v1 + 9 v2 + 16 v12

julia> GeometricAlgebra.show_multivector(stdout, a; parseable=true)
Multivector{2, 0:2}([1, 4, 9, 16])

```
"""
function show_multivector(io::IO, @nospecialize(a);
                          inline=display_options.inline,
                          groupgrades=something(display_options.groupgrades, !ishomogeneous(a)),
                          showzeros=something(display_options.showzeros, ishomogeneous(a) && dimension(a) < 8),
                          eps=display_options.eps,
                          indent=0,
                          parseable=display_options.parseable,
                          compact=display_options.compact,
                          basis_display_style=get_basis_display_style(signature(a)))
	if parseable
		@static if VERSION ≥ v"1.7"
			type = Base.typeinfo_implicit(typeof(a.comps)) ? constructor(a) : typeof(a)
		else
			type = typeof(a)
		end
		print(io, type, "(", a.comps, ")")
	else
		if groupgrades
			if issmall(a, eps)
				print(io, " "^indent, numberzero(isconcretetype(eltype(a)) ? eltype(a) : Int))
				return
			end

			firstgroup = true
			for k in grade(a)
				ak = grade(a, k)

				if !showzeros
					issmall(ak, eps) && continue
				end
				firstgroup || print(io, inline ? " + " : "\n")
				if firstgroup || !inline
					print(io, " "^indent)
				end

				inline && print(io, "(")
				show_multivector_row(io, ak; showzeros, compact, eps, basis_display_style)
				inline && print(io, ")")

				firstgroup = false
			end
		else
			if inline
				show_multivector_row(io, a; indent, showzeros, compact, eps, basis_display_style)
			else
				show_multivector_col(io, a; indent, showzeros, compact, eps, basis_display_style)
			end
		end
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


function Base.show(io::IO, @nospecialize(a::BasisBlade))
	if :__PRETTY_TABLES_DATA__ ∈ keys(io)
		# printing within human-readable table
		show_blade(IOContext(io, :color=>false, :compact=>true), a)
	else
		show_blade(io, a, parseable=true)
	end
end
function Base.show(io::IO, ::MIME"text/plain", @nospecialize(a::BasisBlade))
	if :typeinfo ∈ keys(io) 
		show_blade(io, a)
	else
		show_header(io, a)
		print(io, " ")
		show_blade(io, a)
	end
end

function Base.show(io::IO, @nospecialize(a::Multivector))
	if :__PRETTY_TABLES_DATA__ ∈ keys(io)
		# printing within human-readable table
		show_multivector(IOContext(io, :color=>false), a; inline=true, compact=true, showzeros=false)
	else
		show_multivector(io, a, parseable=true)
	end
end
function Base.show(io::IO, ::MIME"text/plain", @nospecialize(a::Multivector))
	if :typeinfo ∈ keys(io)
		show_multivector(io, a; inline=true, compact=true, showzeros=false)
	else
		show_header(io, a)
		show_multivector(io, a; indent=1)
	end
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
