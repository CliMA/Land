# function to get leaf p, tree dedp, leaf a, and leaf ci for single layered tree, test only
function Yujie111GetPACGT(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            flow,
            leaf,
            envir::AirLayer{FT},
            displaying=false
) where {FT<:AbstractFloat}
    # 0. unpack the data
    r_sl, q_sl, e_sl, r_sh, q_sh, e_sh = leaf
    @unpack t_air, p_a, vpd, wind = envir;

    q_all = q_sl*r_sl + q_sh*r_sh
    e_all = e_sl + e_sh

    # 1. get leaf Ps
    p_leaf = Yujie111GetP(node, flow)
    p_crit = -xylem_p_crit(node.vc_leaf);

    # 2. get A for combined layer
    v_air  = saturation_vapor_pressure(t_air + FT(K_0)) / 1000
    e_flow = flow * FT(KG_H_2_MOL_S)
    t_leaf = max(-50, Yujie111GetLeafTem(node, t_air, e_all, e_flow, false, wind));
    leaf_temperature_dependence!(photo_set, node.ps, envir, t_leaf+FT(K_0));
    d_leaf = node.ps.p_sat - envir.p_H₂O;
    if (d_leaf > 0) && (t_leaf>0) && (p_leaf<p_crit)
        g_lw   = e_flow / node.laba / d_leaf * envir.p_atm
        g_lc   = g_lw / 1.6
        if g_lw < node.g_sla * relative_diffusive_coefficient(FT(K_0) + (t_leaf+t_air)/2)
            leaf_photo_from_glc!(photo_set, node.ps, envir, g_lc);
            a_leaf = node.ps.An;
            c_leaf = node.ps.p_i;
        else
            p_leaf = -Inf
            a_leaf = -Inf
            c_leaf = -Inf
            g_lw   = -Inf
            t_leaf = -Inf
        end
    else
        p_leaf = -Inf
        a_leaf = -Inf
        c_leaf = -Inf
        g_lw   = -Inf
        t_leaf = -Inf
    end

    # last step return values
    return [flow p_leaf a_leaf c_leaf g_lw t_leaf]
end
