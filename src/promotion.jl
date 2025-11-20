#= Promotion =#

for op in [
	:(Base.:+),
	:geometric_prod,
	:scalar_prod,
]
	@eval $op(a::AbstractMultivector...) = $op(promote(a...)...)
end
graded_prod(fn, a::AbstractMultivector...) = graded_prod(fn, promote(a...)...)


"""
	signature_promote_rule(Val(Sig1), Val(Sig2)) -> Sig

Define a promotion rule for multivector signatures.

This works in much the same way as `Base.promote_rule`, except that this does not
operate on the full types, but only on the metric signature parameter of `AbstractMultivector`.

Operations between multivectors of different signatures is not defined unless there is an applicable
method of [`signature_promote_rule`](@ref) defining the common signature to convert arguments to.
Signatures must be wrapped in `Val{}`, not `Type{}`, because signatures are allowed to be `isbits`
values as well as types.

To convert multivectors between signatures, there must be a corresponding method `signature_convert(::Val{Sig}, a)`.

See also [`signature_convert`](@ref).

# Example

```julia
signature_promote_rule(::Val{CGA{Sig}}, ::Val{Sig}) where Sig = CGA{Sig}
```
"""
function signature_promote_rule end
signature_promote_rule(::Val, ::Val) = nothing

"""
	signature_convert(::Val{Sig}, a::AbstractMultivector{Sig′}) -> AbstractMultivector{Sig}

Convert a multivector to a different metric signature.

If it makes sense for multivectors of specific distinct signatures to be interoperable,
this method should be defined along with a signature promotion rule.

See also [`signature_promote_rule`](@ref).


# Example

```julia
signature_convert(::Val{CGA{Sig}}, a::AbstractMultivector{Sig}) where Sig = embed(CGA{Sig}, a)
```
"""
signature_convert(::Val{Sig}, a::AbstractMultivector{Sig}) where {Sig} = a

signature_convert(::Val{Sig}, a::AbstractMultivector{Sig′}) where {Sig,Sig′} = error("""
Cannot convert multivector with signature $Sig′ to signature $Sig.

Hint: Define the method `$signature_convert(::Val{$Sig}, a::AbstractMultivector{$Sig′})` to enable promotion.
""")

function Base.promote_rule(T::Type{<:AbstractMultivector{Sig}}, T′::Type{<:AbstractMultivector{Sig′}}) where {Sig,Sig′}
	sig = @something(
		signature_promote_rule(Val(Sig), Val(Sig′)),
		signature_promote_rule(Val(Sig′), Val(Sig)),
		Some(nothing)
	)
	isnothing(sig) && return Union{}
	AbstractMultivector{sig}
end

function promote_signature(::Val{Sig}, ::Val{Sig′}) where {Sig,Sig′}
	signature(promote_type(AbstractMultivector{Sig}, AbstractMultivector{Sig′}))
end

function Base.convert(::Type{AbstractMultivector{Sig}}, a::AbstractMultivector{Sig′}) where {Sig,Sig′}
	signature_convert(Val(Sig), a)
end

# this is the error shown when promotion fails
function Base.sametype_error(x::Tuple{Vararg{AbstractMultivector}})
	error("""
	Could not promote multivectors of signature $(join(signature.(x), ", ", " and ")).

	Hint: If appropriate, define methods for `signature_promote_rule()` and
	`signature_convert()` to enable promotion.
	""")
end
