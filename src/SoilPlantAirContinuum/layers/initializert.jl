"""
    initialize_spac_canopy!(node::SPACMono{FT}) where {FT<:AbstractFloat}

Initialize the RT parameters for a given
- `node` [`SPACMono`](@ref) type struct
"""
function initialize_spac_canopy!(
            node::SPACMono{FT}
) where {FT<:AbstractFloat}
    # 0.1 create variables required
    @unpack n_canopy, wl_set = node;
    canopy_rt = node.canopy_rt;
    fraction_sl::Array{FT,1} = repeat(canopy_rt.lidf, outer=[ canopy_rt.nAzi ]) / length(canopy_rt.lazitab);
    n_sl = length(canopy_rt.lidf) * length(canopy_rt.lazitab);

    # fluspect the canopy layers
    for ican in 1:n_canopy
        fluspect!(node.leaves_rt[ican], wl_set);
    end

    # Four Different steps to compute Short-Wave RT
    canopy_geometry!(canopy_rt, node.angles, node.can_opt, node.rt_con)
    canopy_matrices!(node.leaves_rt, node.can_opt);
    short_wave!(canopy_rt, node.can_opt, node.can_rad, node.in_rad, node.soil_opt, node.rt_con);
    canopy_fluxes!(canopy_rt, node.can_opt, node.can_rad, node.in_rad, node.soil_opt, node.leaves_rt, wl_set, node.rt_con);

    # Compute Long Wave (Last term is LW incoming in W m^-2)
    thermal_fluxes!(node.leaves_rt, node.can_opt, node.can_rad, canopy_rt, node.soil_opt, [FT(400.0)], wl_set);

    # update the canopy leaf area partition information
    for i_can in 1:n_canopy
        rt_layer = n_canopy + 1 - i_can;

        # calculate the fraction of sunlit and shaded leaves
        f_view = (node.can_opt.Ps[rt_layer] + node.can_opt.Ps[rt_layer+1]) / 2;

        for i_leaf in 1:n_sl
            node.plant_ps[i_can].APAR[i_leaf]  = node.can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
            node.plant_ps[i_can].APAR[i_leaf] += node.can_rad.absPAR_sunCab[(rt_layer-1)*n_sl + i_leaf] * FT(1e6);
            node.plant_ps[i_can].LAIx[i_leaf]  = f_view * fraction_sl[i_leaf];
        end
        node.plant_ps[i_can].APAR[end] = node.can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
        node.plant_ps[i_can].LAIx[end] = 1 - f_view;

        update_leaf_TP!(node.photo_set, node.plant_ps[i_can], node.plant_hs.leaves[i_can], node.envirs[i_can]);

        node.envirs[i_can].t_air = K_25(FT);
        node.envirs[i_can].p_sat = saturation_vapor_pressure(K_25(FT));
        node.envirs[i_can].p_H₂O = node.envirs[i_can].p_sat / 2;
    end

    return nothing
end
