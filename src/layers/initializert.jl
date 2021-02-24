###############################################################################
#
# Initialize SPACMono structure
#
###############################################################################
"""
    initialize_spac_canopy!(node::SPACMono{FT}) where {FT<:AbstractFloat}

Initialize the RT parameters for a given
- `node` [`SPACMono`](@ref) type struct
"""
function initialize_spac_canopy!(
            node::SPACMono{FT}
) where {FT<:AbstractFloat}
    # 0.1 create variables required
    @unpack angles, can_opt, can_rad, canopy_rt, envirs, in_rad, leaves_rt,
            n_canopy, plant_hs, plant_ps, photo_set, rt_con, soil_opt,
            wl_set = node;
    fraction_sl::Array{FT,1} = repeat(canopy_rt.lidf, outer=[canopy_rt.nAzi]) /
                               length(canopy_rt.lazitab);
    n_sl = length(canopy_rt.lidf) * length(canopy_rt.lazitab);

    # fluspect the canopy layers
    for ican in 1:n_canopy
        fluspect!(leaves_rt[ican], wl_set);
    end

    # Four Different steps to compute Short-Wave RT
    canopy_geometry!(canopy_rt, angles, can_opt, rt_con)
    canopy_matrices!(leaves_rt, can_opt);
    short_wave!(canopy_rt, can_opt, can_rad, in_rad,
                soil_opt, rt_con);
    canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad,
                   soil_opt, leaves_rt, wl_set, rt_con);

    # Compute Long Wave (Last term is LW incoming in W m^-2)
    thermal_fluxes!(leaves_rt, can_opt, can_rad, canopy_rt,
                    soil_opt, [FT(400.0)], wl_set);

    # update the canopy leaf area partition information
    for i_can in 1:n_canopy
        rt_layer = n_canopy + 1 - i_can;
        iPS      = plant_ps[i_can];

        # calculate the fraction of sunlit and shaded leaves
        f_view = (can_opt.Ps[rt_layer] + can_opt.Ps[rt_layer+1]) / 2;

        for iLF in 1:n_sl
            iPS.APAR[iLF] = can_rad.absPAR_sunCab[(rt_layer-1)*n_sl + iLF] *
                            FT(1e6);
            iPS.LAIx[iLF] = f_view * fraction_sl[iLF];
        end
        iPS.APAR[end] = can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
        iPS.LAIx[end] = 1 - f_view;

        update_leaf_TP!(photo_set, iPS, plant_hs.leaves[i_can], envirs[i_can]);

        envirs[i_can].t_air = K_25(FT);
        envirs[i_can].p_sat = saturation_vapor_pressure(K_25(FT));
        envirs[i_can].p_H₂O = envirs[i_can].p_sat / 2;
    end

    return nothing
end
