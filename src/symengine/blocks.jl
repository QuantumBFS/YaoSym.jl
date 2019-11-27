using YaoBlocks, SymEngine, LuxurySparse, LinearAlgebra

Base.promote_rule(::Type{Bool}, ::Type{Basic}) = Basic
Base.conj(x::Basic) = real(x) - imag(x)

const SymReal = Union{Basic, SymEngine.BasicRealNumber}
YaoBlocks.RotationGate(block::GT, theta::T) where {N, T <: SymReal, GT<:AbstractBlock{N}} = RotationGate{N, T, GT}(block, theta)

YaoBlocks.phase(θ::SymReal) = PhaseGate(θ)
YaoBlocks.shift(θ::SymReal) = ShiftGate(θ)

YaoBlocks.mat(::Type{Basic}, ::HGate) = 1/sqrt(Basic(2)) * Basic[1 1;1 -1]
YaoBlocks.mat(::Type{Basic}, ::XGate) = Basic[0 1;1 0]
YaoBlocks.mat(::Type{Basic}, ::YGate) = Basic[0 -1im;1im 0]
YaoBlocks.mat(::Type{Basic}, ::ZGate) = Basic[1 0;0 -1]

YaoBlocks.mat(gate::ShiftGate{<:SymReal}) =
    Diagonal([1.0, exp(im * gate.theta)])
YaoBlocks.mat(gate::PhaseGate{<:SymReal}) =
    exp(im * gate.theta) * IMatrix{2}()
function YaoBlocks.mat(R::RotationGate{N, <:SymReal}) where N
    I = IMatrix{1<<N}()
    return I * cos(R.theta / 2) - im * sin(R.theta / 2) * mat(Basic,R.block)
end

YaoBlocks.PSwap{N}(locs::Tuple{Int, Int}, θ::SymReal) where N = YaoBlocks.PutBlock{N}(rot(ConstGate.SWAPGate(), θ), locs)

YaoBlocks.pswap(n::Int, i::Int, j::Int, α::SymReal) = PSwap{n}((i,j), α)
YaoBlocks.pswap(i::Int, j::Int, α::SymReal) = n->pswap(n,i,j,α)

export subs, chiparam
chiparam(blk::RotationGate, param) = rot(blk.block, param)
chiparam(blk::ShiftGate, param) = shift(param)
chiparam(blk::PhaseGate, param) = phase(param)
chiparam(blk::TimeEvolution, param) = time_evolve(blk.H, param, tol=blk.tol)
chiparam(blk::AbstractBlock, params...) = niparams(blk) == length(params) == 0 ? blk : NotImplementedError(chiparam, blk, params...)
SymEngine.subs(blk::AbstractBlock, args...; kwargs...) = subs(Basic, blk, args...; kwargs...)
SymEngine.subs(::Type{T}, blk::RotationGate, args...; kwargs...) where T = rot(blk.block, T(subs(blk.theta, args...; kwargs...)))
SymEngine.subs(::Type{T}, blk::ShiftGate, args...; kwargs...) where T = rot(blk.block, T(subs(blk.theta, args...; kwargs...)))
SymEngine.subs(::Type{T}, blk::ShiftGate, args...; kwargs...) where T = rot(blk.block, T(subs(blk.theta, args...; kwargs...)))
SymEngine.subs(::Type{T}, blk::TimeEvolution, args...; kwargs...) where T = time_evolve(blk.H, T(subs(blk.dt, args...; kwargs...)), tol=blk.tol)
function SymEngine.subs(c::AbstractBlock, args...; kwargs...)
    nparameters(c) == 0 && return c
    if niparams(c) == 0
        chsubblocks(c, [subs(blk, args..., kwargs...) for blk in subblocks(blk)])
    else
        throw(NotImplementedError(:subs, (c, args...)))
    end
end
