# function to get leaf p, tree dedp, leaf a, and leaf ci
function Yujie111GetPACGTs(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            f_sl,
            f_sh,
            displaying=false
) where {FT<:AbstractFloat}
    # unpack the data
    @unpack frac_sh, frac_sl, la_sh, la_sl, par_sh, par_sl, rad_sh, rad_sl = node.container2L;

    # calculate gas exchangr for sunlit and shaded layers
    leaf_gas_exchange!(node, photo_set, f_sl, par_sl, rad_sl, la_sl, node.container2L.cont_sl);
    leaf_gas_exchange!(node, photo_set, f_sh, par_sh, rad_sh, la_sh, node.container2L.cont_sh);

    return nothing
end
