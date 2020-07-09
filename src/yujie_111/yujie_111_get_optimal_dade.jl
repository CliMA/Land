function Yujie111GetOptimaldAdE(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            zenith=30.0,
            r_all=1000.0
) where {FT<:AbstractFloat}
    # partition the sunlit and shade layer
    canopy    = Yujie111GetLeafPartition(node, zenith, r_all)

    # get the optimal for this partition
    f_sl,f_sh = Yujie111GetOptimalFs(node, photo_set, envir, zenith, r_all)
    tmp_re    = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
    anet      = tmp_re[1,3]*canopy[1] + tmp_re[2,3]*canopy[4]

    # calculate dade
    de_sl     = f_sl/(f_sl+f_sh)
    de_sh     = f_sh/(f_sl+f_sh)
    f_sl     += de_sl
    f_sh     += de_sh
    tmp_de    = Yujie111GetPACGTs(node, photo_set, f_sl, f_sh, canopy, envir)
    anet_de   = tmp_de[1,3]*canopy[1] + tmp_de[2,3]*canopy[4]
    dade      = anet_de - anet
    return dade
end