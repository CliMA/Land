# ∂A∂E = ∂Θ∂E, use ∂Θ∂E to calculate this?




function Yujie111GetOptimaldAdE(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            zenith=30.0,
            r_all=1000.0
) where {FT<:AbstractFloat}
    # partition the sunlit and shade layer
    big_leaf_partition!(node, zenith, r_all)
    @unpack frac_sh, frac_sl = node.container2L;

    # get the optimal for this partition
    Yujie111GetOptimalFs(node, photo_set, zenith, r_all)
    Yujie111GetPACGTs(node, photo_set, node.opt_f_sl, node.opt_f_sh)
    anet    = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;

    # calculate dade
    de_sl   = node.opt_f_sl/(node.opt_f_sl+node.opt_f_sh)
    de_sh   = 1 - de_sl
    f_sl    = node.opt_f_sl + de_sl * FT(0.01)
    f_sh    = node.opt_f_sl + de_sh * FT(0.01)
    Yujie111GetPACGTs(node, photo_set, f_sl, f_sh)
    anet_de = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
    dade    = (anet_de - anet) * 100
    return dade
end
