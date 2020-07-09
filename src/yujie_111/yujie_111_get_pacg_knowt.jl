# function to test birch dataset, know t_leaf
function Yujie111GetPACGKnowT(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            conditions,
            flow
) where {FT<:AbstractFloat}
    # 0. unpack data




    par,p_co2,vp_air,t_leaf = conditions





    # 1. get leaf p
    p_leaf = Yujie111GetP(node, flow)

    # 2. get A for combined layer
    e_flow = flow * FT(KG_H_2_MOL_S)
    leaf_temperature_dependence!(photo_set, node.ps, envir, t_leaf+FT(K_0));
    d_leaf = node.ps.p_sat - envir.p_H₂O
    g_lw   = e_flow / node.laba / d_leaf * envir.p_atm
    g_lc   = g_lw / 1.6
    leaf_photo_from_glc!(photo_set, node.ps, envir, g_lc);

    a_leaf = node.ps.An;
    c_leaf = node.ps.p_i;

    # last step return values
    return [p_leaf a_leaf c_leaf g_lw]
end
