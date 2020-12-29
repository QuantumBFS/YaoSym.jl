using ..SymEngine
import YaoBase: rot_mat

rot_mat(::Type{T}, ::Val{:Rx}, theta::Basic) where {T} =
    Basic[cos(theta / 2) -im * sin(theta / 2); -im * sin(theta / 2) cos(theta / 2)]
rot_mat(::Type{T}, ::Val{:Ry}, theta::Basic) where {T} =
    Basic[cos(theta / 2) -sin(theta / 2); sin(theta / 2) cos(theta / 2)]
rot_mat(::Type{T}, ::Val{:Rz}, theta::Basic) where {T} =
    Diagonal(Basic[exp(-im * theta / 2), exp(im * theta / 2)])
rot_mat(::Type{T}, ::Val{:CPHASE}, theta::Basic) where {T} = Diagonal(Basic[1, 1, 1, exp(im * theta)])
rot_mat(::Type{T}, ::Val{:PSWAP}, theta::Basic) where {T} = rot_mat(Basic, Const.SWAP, theta)

for G in [:(Val{:Rx}), :(Val{:Ry}), :(Val{:Rz}), :(Val{:CPHASE}), :(Val)]
    # forward single gates
    @eval function YaoBase.instruct!(
        state::AbstractVecOrMat{T},
        g::$G,
        locs::Tuple{Int},
        control_locs::NTuple{N1,Int},
        control_bits::NTuple{N2,Int},
        theta::Basic,
    ) where {T,N1,N2}
        m = rot_mat(T, g, theta)
        instruct!(state, m, locs, control_locs, control_bits)
    end
    # forward single gates
    @eval function YaoBase.instruct!(
        state::AbstractVecOrMat{T},
        g::$G,
        locs::Tuple{Int},
        theta::Basic,
    ) where {T,N1}
        instruct!(state, g, locs, (), (), theta)
    end
end

@eval function YaoBase.instruct!(
    state::AbstractVecOrMat{T},
    g::Val{:PSWAP},
    locs::Tuple{Int,Int},
    theta::Basic,
) where {T,N1}
    instruct!(state, g, locs, (), (), theta)
end
