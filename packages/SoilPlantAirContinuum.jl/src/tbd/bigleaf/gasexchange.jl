#=
###############################################################################
#
# Calculate gas exchange from flow rate only for required parameters
#
###############################################################################
"""
    leaf_gas_exchange_nonopt!(
                node::SPACSimple{FT},
                flow::FT,
                par::FT,
                rad::FT,
                la::FT,
                container::SPACContainer1L{FT}
    ) where {FT<:AbstractFloat}
    leaf_gas_exchange_nonopt!(
                node::SPACSimple{FT},
                flow::FT
    ) where {FT<:AbstractFloat}
    leaf_gas_exchange_nonopt!(
                node::SPACSimple{FT},
                f_sl::FT,
                f_sh::FT
    ) where {FT<:AbstractFloat}

Simulate leaf level gas exchange and fill it into the `container` for 1-layer
    or 2-layer canopy, given
- `node` [`SPACSimple`] type struct
- `flow` Flow rate per basal area into the leaves (e.g., for sunlit leaves)
- `f_sl` Flow rate per basal area into the sunlit leaves
- `f_sh` Flow rate per basal area into the shaded leaves
- `par` Leaf-level photosynthetic active radiation
- `rad` Leaf-level absorbed radiative energy
- `la` Leaf area of the leaves (total or each layer)
- `container` [`SPACContainer1L`] type container
"""
function leaf_gas_exchange_nonopt!(
            node::SPACSimple{FT},
            flow::FT,
            par::FT,
            rad::FT,
            la::FT,
            container::SPACContainer1L{FT}
) where {FT<:AbstractFloat}
    # if flow < 0
    if flow < 0
        container.an = FT(-Inf);

    # if flow >= 0
    else
        # 0. unpack required variables
        @unpack envir, g_max, width = node;
        @unpack P_AIR, p_H₂O, wind = envir;

        # 1. calculate leaf temperature from the flow rate
        t_leaf = max(200, leaf_temperature(node, rad, flow));

        # 2. update leaf photosynthetic variables and leaf-to-air VPD
        node.ps.apar = par;
        photosystem_temperature_dependence!(node.ps.PSM, envir, t_leaf);
        node.ps.p_H₂O_sat = saturation_vapor_pressure(t_leaf);
        d_leaf = node.ps.p_H₂O_sat - p_H₂O;

        # 3. update f_vis and f_st in leaf and calculate water potentials
        # TODO do not account for temperature effects for flow now
        # More reasonable functions need to be added
        # vc_temperature_effects!(node.hs.leaf, t_leaf);

        # if flow == 0
        if flow == 0
            g_lw = FT(0);
            g_lc = FT(1e-6);
            leaf_photosynthesis!(node.ps, envir, GCO₂Mode(), g_lc);
            container.an = node.ps.PSM.a_net;

        # if flow > 0 and reasonable
        elseif (d_leaf > 0) && (t_leaf > T_0(FT))
            t_cor = relative_diffusive_coefficient( (t_leaf+node.envir.t)/2 );
            g_lw  = flow / la / d_leaf * P_AIR;
            g_bw  = boundary_layer_conductance(wind, width);
            g_bc  = g_bw / FT(1.35);
            g_sw  = 1 / (1/g_lw - 1/(g_bw*t_cor));
            g_sc  = g_sw / FT(1.6);
            g_lc  = max(FT(1e-6), 1 / (1/g_bc + 1/g_sc));
            g_lim = 1 / (1/g_bw + 1/g_max);
            if g_lw < g_lim * t_cor
                leaf_photosynthesis!(node.ps, envir, GCO₂Mode(), g_lc);
                container.an = node.ps.PSM.a_net;
            else
                container.an = FT(-Inf);
            end

        # if flow > 0 and unrealistic
        else
            container.an = FT(-Inf);
        end
    end

    return nothing
end




function leaf_gas_exchange_nonopt!(
            node::SPACSimple{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # unpack the data
    @unpack frac_sh, frac_sl, par_sh, par_sl, rad_sh, rad_sl = node.container2L;

    # calculate mean par and rad per leaf area, then gas exchange rate
    par_mean = par_sl * frac_sl + par_sh * frac_sh;
    rad_mean = rad_sl * frac_sl + rad_sh * frac_sh;
    leaf_gas_exchange_nonopt!(node, flow, par_mean, rad_mean, node.laba, node.container1L);

    node.containerOP = (node.ec - flow) * (node.container1L).an;

    return nothing
end




function leaf_gas_exchange_nonopt!(
            node::SPACSimple{FT},
            f_sl::FT,
            f_sh::FT
) where {FT<:AbstractFloat}
    # unpack the data
    @unpack frac_sh, frac_sl, la_sh, la_sl, par_sh, par_sl, rad_sh,
            rad_sl = node.container2L;

    # calculate gas exchangr for sunlit and shaded layers
    leaf_gas_exchange_nonopt!(node, f_sl, par_sl, rad_sl, la_sl, node.container2L.cont_sl);
    leaf_gas_exchange_nonopt!(node, f_sh, par_sh, rad_sh, la_sh, node.container2L.cont_sh);

    a_sum = frac_sl * node.container2L.cont_sl.an + frac_sh * node.container2L.cont_sh.an;
    e_sum = f_sl + f_sh;
    node.containerOP = (node.ec - e_sum) * a_sum;

    return nothing
end








###############################################################################
#
# Calculate gas exchange from flow rate considering temperature dynamics
#
###############################################################################
"""
    leaf_gas_exchange!(
                node::SPACSimple{FT},
                flow::FT,
                par::FT,
                rad::FT,
                la::FT,
                container::SPACContainer1L{FT}
    ) where {FT<:AbstractFloat}
    leaf_gas_exchange!(
                node::SPACSimple{FT},
                flow::FT
    ) where {FT<:AbstractFloat}
    leaf_gas_exchange!(
                node::SPACSimple{FT},
                f_sl::FT,
                f_sh::FT
    ) where {FT<:AbstractFloat}

Simulate leaf level gas exchange and fill it into the `container` for 1-layer
    or 2-layer canopy, given
- `node` [`SPACSimple`] type struct
- `flow` Flow rate per basal area into the leaves (e.g., for sunlit leaves)
- `f_sl` Flow rate per basal area into the sunlit leaves
- `f_sh` Flow rate per basal area into the shaded leaves
- `par` Leaf-level photosynthetic active radiation
- `rad` Leaf-level absorbed radiative energy
- `la` Leaf area of the leaves (total or each layer)
- `container` [`SPACContainer1L`] type container
"""
function leaf_gas_exchange!(
            node::SPACSimple{FT},
            flow::FT,
            par::FT,
            rad::FT,
            la::FT,
            container::SPACContainer1L{FT}
) where {FT<:AbstractFloat}
    # 0. unpack required variables
    @unpack envir = node;
    @unpack P_AIR, p_H₂O = envir;

    # 1. calculate leaf temperature from the flow rate
    t_leaf = max(200, leaf_temperature(node, rad, flow));

    # 2. update leaf photosynthetic variables and leaf-to-air VPD
    node.ps.apar = par;
    photosystem_temperature_dependence!(node.ps.PSM, envir, t_leaf);
    node.ps.p_H₂O_sat = saturation_vapor_pressure(t_leaf);
    d_leaf = node.ps.p_H₂O_sat - p_H₂O;

    # 3. update f_vis and f_st in leaf and calculate water potentials
    # TODO do not account for temperature effects for flow now
    # More reasonable functions need to be added
    # vc_temperature_effects!(node.hs.leaf, t_leaf);
    p_leaf = end_pressure(node.hs, flow);

    # 4. calculate photosynthesis
    g_lw = flow / la / d_leaf * P_AIR;
    g_lc = max(FT(1e-6), g_lw / FT(1.6));
    leaf_photosynthesis!(node.ps, envir, GCO₂Mode(), g_lc);
    container.ag = node.ps.PSM.a_gross;
    container.an = node.ps.PSM.a_net;
    container.c  = node.ps.p_CO₂_i;
    container.e  = flow;
    container.gh = g_lw;
    container.p  = p_leaf;
    container.t  = t_leaf;

    return nothing
end




function leaf_gas_exchange!(
            node::SPACSimple{FT},
            flow::FT
) where {FT<:AbstractFloat}
    # unpack the data
    @unpack container2L, laba = node;
    @unpack frac_sh, frac_sl, par_sh, par_sl, rad_sh, rad_sl = container2L;

    # calculate mean par and rad per leaf area, then gas exchange rate
    par_mean = par_sl * frac_sl + par_sh * frac_sh;
    rad_mean = rad_sl * frac_sl + rad_sh * frac_sh;
    leaf_gas_exchange!(node, flow, par_mean, rad_mean, laba, node.container1L);

    return nothing
end




function leaf_gas_exchange!(
            node::SPACSimple{FT},
            f_sl::FT,
            f_sh::FT
) where {FT<:AbstractFloat}
    # unpack the data
    @unpack frac_sh, frac_sl, la_sh, la_sl, par_sh, par_sl, rad_sh,
            rad_sl = node.container2L;

    # calculate gas exchangr for sunlit and shaded layers
    leaf_gas_exchange!(node, f_sl, par_sl, rad_sl, la_sl, node.container2L.cont_sl);
    leaf_gas_exchange!(node, f_sh, par_sh, rad_sh, la_sh, node.container2L.cont_sh);

    return nothing
end
=#
