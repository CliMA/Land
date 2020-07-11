# function to get leaf p, tree dedp, leaf a, and leaf ci for single layered tree, test only
function Yujie111GetPACGT(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            flow,
            displaying=false
) where {FT<:AbstractFloat}
    # unpack the data
    @unpack frac_sh, frac_sl, par_sh, par_sl, rad_sh, rad_sl = node.container2L;

    # calculate mean par and rad per leaf area, then gas exchange rate
    par_mean = par_sl * frac_sl + par_sh * frac_sh;
    rad_mean = rad_sl * frac_sl + rad_sh * frac_sh;
    leaf_gas_exchange!(node, photo_set, flow, par_mean, rad_mean, node.laba, node.container1L);

    return nothing
end
