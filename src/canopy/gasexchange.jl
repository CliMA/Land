###############################################################################
#
# Calculate gas exchange from flow rate considering temperature dynamics
#
###############################################################################

function leaf_gas_exchange!(
            node::SPACSimple{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            flow::FT,
            par::FT,
            rad::FT,
            la::FT,
            container::SPACContainer1L{FT}
) where {FT<:AbstractFloat}
    # 0. unpack required variables
    @unpack p_atm, p_H₂O, t_air = node.envir;

    # 1. calculate leaf temperature from the flow rate
    t_leaf = max(200, leaf_temperature(node, rad, flow));

    # 2. update leaf photosynthetic variables and leaf-to-air VPD
    node.ps.APAR = par;
    leaf_temperature_dependence!(photo_set, node.ps, node.envir, t_leaf);
    d_leaf = node.ps.p_sat - p_H₂O;

    # 3. update f_vis and f_st in leaf and calculate water potentials
    vc_temperature_effects!(node.hs.leaf, t_leaf);
    p_leaf = xylem_p_from_flow(node.hs, flow)
    p_crit = node.hs.leaf.p_crt;

    #@show t_leaf;
    #@show d_leaf;
    #@show p_leaf;

    # 4. determine whether to compute photosynthesis
    if flow == 0
        g_lw = FT(0);
        g_lc = FT(1e-6);
        leaf_photo_from_glc!(photo_set, node.ps, node.envir, g_lc);
        container.ag = node.ps.Ag;
        container.an = node.ps.An;
        container.c  = node.ps.p_i;
        container.e  = flow;
        container.gh = g_lw;
        container.p  = p_leaf;
        container.t  = t_leaf;
    elseif (d_leaf > 0) && (t_leaf > K_0) && (p_leaf > p_crit)
        g_lw = flow / la / d_leaf * p_atm
        g_lc = max(1e-6, g_lw / 1.6)
        if g_lw < node.g_sla * relative_diffusive_coefficient( (t_leaf+t_air)/2 )
            leaf_photo_from_glc!(photo_set, node.ps, node.envir, g_lc);
            container.ag = node.ps.Ag;
            container.an = node.ps.An;
            container.c  = node.ps.p_i;
            container.e  = flow;
            container.gh = g_lw;
            container.p  = p_leaf;
            container.t  = t_leaf;
        else
            container.ag = FT(-Inf);
            container.an = FT(-Inf);
            container.c  = FT(-Inf);
            container.e  = FT(-Inf);
            container.gh = FT(-Inf);
            container.p  = FT(-Inf);
            container.t  = FT(-Inf);
        end
    else
        container.ag = FT(-Inf);
        container.an = FT(-Inf);
        container.c  = FT(-Inf);
        container.e  = FT(-Inf);
        container.gh = FT(-Inf);
        container.p  = FT(-Inf);
        container.t  = FT(-Inf);
    end

    return nothing
end