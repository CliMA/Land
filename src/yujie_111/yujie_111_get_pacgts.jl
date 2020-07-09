# function to get leaf p, tree dedp, leaf a, and leaf ci
function Yujie111GetPACGTs(
            node::Yujie111{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            f_sl,
            f_sh,
            leaf,
            envir::AirLayer{FT},
            displaying=false
) where {FT<:AbstractFloat}
    # 0. unpack the data
    r_sl, q_sl, e_sl, r_sh, q_sh, e_sh = leaf
    @unpack t_air, p_a, vpd, wind = envir;

    # 1. get leaf Ps
    p_l_sl, p_l_sh = Yujie111GetPs(node, f_sl, f_sh, r_sl)
    p_crit         = -xylem_p_crit(node.vc_leaf);

    # 2. get A for sunlit layer
    e_l_sl = f_sl * FT(KG_H_2_MOL_S)
    t_l_sl = max(-50, Yujie111GetLeafTem(node, t_air, e_sl, e_l_sl, false, wind));
    leaf_temperature_dependence!(photo_set, node.ps, envir, t_l_sl+FT(K_0));
    d_l_sl = node.ps.p_sat - envir.p_H₂O
    if f_sl==0
        g_h_sl = 0.0
        g_c_sl = 0.0
        leaf_photo_from_glc!(photo_set, node.ps, envir, g_c_sl);
        a_l_sl = node.ps.An;
        c_l_sl = node.ps.p_i;
    elseif (d_l_sl > 0) && (t_l_sl>0) && (p_l_sl<p_crit)
        g_h_sl = e_l_sl / (node.laba*r_sl) / d_l_sl * envir.p_atm
        g_c_sl = g_h_sl / 1.6
        if g_h_sl < node.g_sla * relative_diffusive_coefficient(FT(K_0) + (t_l_sl+t_air)/2)
            leaf_photo_from_glc!(photo_set, node.ps, envir, g_c_sl);
            a_l_sl = node.ps.An;
            c_l_sl = node.ps.p_i;
        else
            p_l_sl = -Inf
            a_l_sl = -Inf
            c_l_sl = -Inf
            g_h_sl = -Inf
            t_l_sl = -Inf
        end
    else
        p_l_sl = -Inf
        a_l_sl = -Inf
        c_l_sl = -Inf
        g_h_sl = -Inf
        t_l_sl = -Inf
    end

    # 3. get A for shade layer
    e_l_sh = f_sh * FT(KG_H_2_MOL_S)
    t_l_sh = max(-50, Yujie111GetLeafTem(node, t_air, e_sh, e_l_sh, true, wind));
    leaf_temperature_dependence!(photo_set, node.ps, envir, t_l_sh+FT(K_0));
    d_l_sh = node.ps.p_sat - envir.p_H₂O
    if f_sh==0
        g_h_sh = 0.0
        g_c_sh = 0.0
        leaf_photo_from_glc!(photo_set, node.ps, envir, g_c_sh);
        a_l_sh = node.ps.An;
        c_l_sh = node.ps.p_i;
    elseif (d_l_sh > 0) && (t_l_sh>0) && (p_l_sh<p_crit)
        g_h_sh = e_l_sh / (node.laba*r_sh) / d_l_sh * envir.p_atm
        g_c_sh = g_h_sh / 1.6
        if g_h_sh < node.g_sla * relative_diffusive_coefficient(FT(K_0) + (t_l_sh+t_air)/2)
            leaf_photo_from_glc!(photo_set, node.ps, envir, g_c_sh);
            a_l_sh = node.ps.An;
            c_l_sh = node.ps.p_i;
        else
            p_l_sh = -Inf
            a_l_sh = -Inf
            c_l_sh = -Inf
            g_h_sh = -Inf
            t_l_sh = -Inf
        end
    else
        p_l_sh = -Inf
        a_l_sh = -Inf
        c_l_sh = -Inf
        g_h_sh = -Inf
        t_l_sh = -Inf
    end


    # last step return values
    return [f_sl p_l_sl a_l_sl c_l_sl g_h_sl t_l_sl;
            f_sh p_l_sh a_l_sh c_l_sh g_h_sh t_l_sh]
end
