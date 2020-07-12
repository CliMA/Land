# function to test birch dataset, know t_leaf
# not tested yet




function Yujie111GetPACGKnowT(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            conditions,
            flow
) where {FT<:AbstractFloat}
    # 0. unpack data




    par,p_co2,vp_air,t_leaf = conditions





    # 1. get leaf p
    p_leaf = xylem_p_from_flow(node.hs, flow)

    # 2. get A for combined layer
    node.ps.APAR = par;
    leaf_temperature_dependence!(photo_set, node.ps, t_leaf+FT(K_0));
    d_leaf = node.ps.p_sat - node.envir.p_H₂O
    g_lw   = flow / node.laba / d_leaf * node.envir.p_atm
    g_lc   = g_lw / 1.6
    leaf_photo_from_glc!(photo_set, node.ps, node.envir, g_lc);

    a_leaf = node.ps.An;
    c_leaf = node.ps.p_i;

    # last step return values
    return [p_leaf a_leaf c_leaf g_lw]
end
