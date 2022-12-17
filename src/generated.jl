"""
	use_symbolic_optim(sig) -> Bool

Whether to use symbolic code generation to optimize operations
in algebras of metric signature `sig`.

By default, this is enabled if `dimension(sig) ≤ 8` as a heuristic
(in many dimensions, algebraic expressions may become too unwieldy).
"""
use_symbolic_optim(sig) = dimension(sig) <= 8


function symbolic_components(label, dims...)
	var = SymbolicUtils.Sym{Array{length(dims),Real}}(label)
	indices = Iterators.product(Base.OneTo.(dims)...)
	Any[SymbolicUtils.Term{Real}(getindex, [var, I...]) for I in indices]
end

symbolic_multivector(a::OrType{<:BasisBlade{Sig,K}}, label) where {Sig,K} = symbolic_multivector(Multivector{Sig,K}, label)
symbolic_multivector(A::OrType{<:Multivector{Sig,K}}, label) where {Sig,K} = Multivector{Sig,K}(symbolic_components(label, ncomponents(A)))

toexpr(a, T) = SymbolicUtils.Code.toexpr(a)
function toexpr(a::AbstractMultivector, T)
	compstype = componentstype(signature(a), ncomponents(a), T)
	comps_expr = SymbolicUtils.Code.toexpr(SymbolicUtils.Code.MakeArray(a.comps, compstype, T))

	:( $(constructor(a))($comps_expr) )
end


"""
	symbolic_multivector_eval(f, x::AbstractMultivector...)

Return `f(x...)` using symbolic optimisation.
This is a generated function which first evaluates `f` on symbolic versions of
the multivector arguments `x` and converts the result to unrolled code.
"""
function _symbolic_multivector_eval(f, x...)
	abc = Symbol.('a' .+ (0:length(x) - 1))
	assignments = [:( $a = components(x[$i]) ) for (i, a) in enumerate(abc)]

	x_symb = symbolic_multivector.(x, abc)
	y_symb = f(x_symb...)

	T = numberorany(promote_type(eltype.(x)...))
	result = toexpr(y_symb, T)

	quote
		let $(assignments...)
			$result
		end
	end
end

@generated symbolic_multivector_eval(f, x...) = _symbolic_multivector_eval(f.instance, x...)

function multivector_eval(f, x::AbstractMultivector{Sig}...) where {Sig}
	if use_symbolic_optim(Sig)
		symbolic_multivector_eval(f, x...)
	else
		f(x...)
	end
end


# way to convert a BasisBlade to a Multivector without allocating a full components array
# TODO: take this more seriously
function components(a::BasisBlade{Sig,K}) where {Sig,K}
	i = findfirst(==(bitsof(a)), componentbits(Val(dimension(Sig)), Val(K)))
	SingletonVector(a.coeff, i, ncomponents(Sig, K))
end
components(a::Multivector) = a.comps

struct SingletonVector{T} <: AbstractVector{T}
	el::T
	index::Int
	length::Int
end
Base.length(a::SingletonVector) = a.length
Base.size(a::SingletonVector) = (length(a),)
Base.getindex(a::SingletonVector{T}, i) where {T} = a.index == i ? a.el : numberzero(T)
function Base.iterate(a::SingletonVector, i = 1)
	i > a.length && return
	(a[i], i + 1)
end
