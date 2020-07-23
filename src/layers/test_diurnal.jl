###############################################################################
#
# Steps of non-steady state stomatal controls
#
###############################################################################
#=
0. Initialize plant and environment with minimal stomatal opening
1. Update transpiration rate
2. Update canopy temperature
3. Update plant hydraulic profile every 10 minutes
4. Calculate assimilation rate
5. Calculate dH₂O/dt and dCO₂/dt for the environment
6. Calculate dGs/dt for the given environment
7. Update stomatal conductance
8. Goto step 1

This function will be splited into many functions later!
=#

function test_diurnal(Δt::FT=FT(10)) where {FT<:AbstractFloat}
    # 0.1 create variables required
    spacmono = SPACMono{FT}();
    @unpack ga, la, n_canopy, n_root = spacmono;
    angles, leaves_rt, can_rt, can_opt, can_rad, soil_opt, in_rad, wl_set = initialize_rt_module(n_layer=n_canopy, LAI=la/ga);
    photo_set   = C3CLM(FT);
    fraction_sl = repeat(can_rt.lidf, outer=[ length(can_rt.lazitab) ]) / length(can_rt.lazitab);
    stoma_mod   = ESMBallBerry{FT}();



    # change in_rad to alter direct and diffuse fraction
    # sum of direct  energy in W:
    # e_dir = sum( (wl_set.swl[2:end] .- wl_set.swl[1:end-1]) .* in_rad.E_direct  ) / 1000
    # sum of diffuse energy in W:
    # e_dif = sum( (wl_set.swl[2:end] .- wl_set.swl[1:end-1]) .* in_rad.E_diffuse ) / 1000
    # target e_dir is tar_dir, do it this way:
    # in_rad.E_direct  .*= tar_dir / e_dir
    # target e_dif is tar_dif, do it this way:
    # in_rad.E_diffuse .*= tar_dif / e_dif



    # 0.2 A diurnal cycle for radiation and Tair
    Tmean  = FT(295.15)
    DeltaT = 3
    omega  = 2 * π / (24*3600)
    t      = range(0, stop=24*3600*2, step=Δt)
    phi_t  = omega * t - π * ones(FT, size(t)) / 2
    PARr_t = zeros(FT, size(t))
    Tair_t = zeros(FT, size(t))
    Tlef_t = zeros(FT, size(t))
    N      = length(PARr_t)

    # 0.3 customize dirnal cycles of temperatures and PAR
    for i_tim in 1:N
        PARr_t[i_tim] = max( sin(phi_t[i_tim]), FT(0.0) )
        Tair_t[i_tim] = Tmean + DeltaT * sin(phi_t[i_tim] - π/3)
        Tlef_t[i_tim] = Tmean + DeltaT * sin(phi_t[i_tim] - π/3)
    end

    # 0.4 initialize the canopy temperature
    for i_can in 1:n_canopy
        # Assume RH = 50%
        leaf_temperature_dependence!(photo_set, spacmono.plant_ps[i_can].ps, spacmono.envirs[i_can], Tlef_t[1]);
        spacmono.envirs[i_can].p_H₂O = spacmono.envirs[i_can].p_sat / 2;
        # latent heat flux as well
        spacmono.plant_ps[i_can].LV = latent_heat_vapor(Tlef_t[i_can]) * 1000 / 18;
        # vc dependecy as well
        vc_temperature_effects!(spacmono.plant_hs.leaves[i_can], Tlef_t[i_can]);
    end

    # run dirnal cycle
    # NPP per ground area
    Agpp_t = zeros(FT, size(t));
    # NPP per ground area
    Asum_t = zeros(FT, size(t));
    # trasnpiration rate per ground area
    Esum_t = zeros(FT, size(t));
    for i_tim in 1:N
        # 0. recalculate the canopy RT stuff
        zenith = zenith_angle(spacmono.latitude, 180 + t[i_tim]/86400 );
        zenith = min(88, zenith);
        angles.tts = zenith;
        canopy_geometry!(can_rt, angles, can_opt);
        canopy_matrices!(leaves_rt, can_opt);
        short_wave!(can_rt, can_opt, can_rad, in_rad, soil_opt);
        canopy_fluxes!(can_rt, can_opt, can_rad, in_rad, soil_opt, leaves_rt, wl_set);
        thermal_fluxes!(leaves_rt, can_opt, can_rad, can_rt, soil_opt, [FT(400.0)], wl_set);
        sif_fluxes!(leaves_rt, can_opt, can_rad, can_rt, soil_opt, wl_set);

        # 2. Update transpiration rate
        a_gpp = FT(0);
        a_sum = FT(0);
        e_sum = FT(0);
        for i_can in 1:n_canopy
            rt_layer = spacmono.n_canopy + 1 - i_can;

            # set the diffuse par to all the leaves
            spacmono.plant_ps[i_can].APAR           .= PARr_t[i_tim]  * can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
            spacmono.plant_ps[i_can].APAR[1:end-1] .+= PARr_t[i_tim] .* reshape(can_rad.absPAR_sunCab[:,:,rt_layer],(:,1))[:,1] .* FT(1e6);

            # calculate the fraction of sunlit and shaded leaves
            f_view = (can_opt.Ps[rt_layer] + can_opt.Ps[rt_layer+1]) / 2;
            la_new = [f_view .* fraction_sl; 1-f_view];
            spacmono.plant_ps[i_can].LAIx .= la_new;

            # update leaf temperature
            # update the parameters from stomatal model
            spacmono.plant_ps[i_can].T = Tlef_t[i_tim];
            update_leaf_TP!(photo_set, spacmono.plant_ps[i_can], spacmono.plant_hs.leaves[i_can], spacmono.envirs[i_can]);
            # uncomment this for Sperry, Eller, and Gentine models
            #update_leaf_AK!(photo_set, spacmono.plant_ps[i_can], spacmono.plant_hs.leaves[i_can], spacmono.envirs[i_can]);

            #canopyi.t_list[1:end-1] .= reshape(can_rad.T_sun3D[:,:,rt_layer],(:,1))[:,1]
            #canopyi.t_list[end]      = can_rad.T_shade[rt_layer]
            # Assume RH = 50%
            spacmono.envirs[i_can].t_air = Tair_t[i_tim];
            spacmono.envirs[i_can].p_sat = saturation_vapor_pressure(Tair_t[i_tim]);
            spacmono.envirs[i_can].p_H₂O = spacmono.envirs[i_can].p_sat / 2;

            # update the flow rates
            g_i_can = sum( spacmono.plant_ps[i_can].Ag .* spacmono.plant_ps[i_can].LAIx ) * spacmono.plant_ps[i_can].LA;
            a_i_can = sum( spacmono.plant_ps[i_can].An .* spacmono.plant_ps[i_can].LAIx ) * spacmono.plant_ps[i_can].LA;
            e_i_can = sum( spacmono.plant_ps[i_can].g_lw .* (spacmono.plant_ps[i_can].p_sat - spacmono.envirs[i_can].p_H₂O) ./ spacmono.envirs[i_can].p_atm .* spacmono.plant_ps[i_can].LAIx ) * spacmono.plant_ps[i_can].LA;

            # calculate the photosynthetic rates
            update_leaf_from_gsw!(photo_set, spacmono.plant_ps[i_can], spacmono.envirs[i_can]);
            leaf_gsw_control!(photo_set, spacmono.plant_ps[i_can], spacmono.envirs[i_can]);

            # use the ball-berry model here for now as the ∂A/∂E and ∂A/∂Θ functions are not yet ready
            gsw_ss = empirical_gsw_from_model(stoma_mod, spacmono.plant_ps[i_can], spacmono.envirs[i_can], FT(1));

            # assume τ = 10 minutes
            spacmono.plant_ps[i_can].g_sw .+= (gsw_ss .- spacmono.plant_ps[i_can].g_sw) / 600 * Δt;

            a_gpp += g_i_can;
            a_sum += a_i_can;
            e_sum += e_i_can;
        end

        Agpp_t[i_tim] = a_gpp / spacmono.ga;
        Asum_t[i_tim] = a_sum / spacmono.ga;
        Esum_t[i_tim] = e_sum / spacmono.ga;

        # update pressure profile every 10 minutes
        # TO BE DONE within PlantHydraulics

        _par = PARr_t[i_tim]  * can_rad.absPAR_shadeCab[1] * FT(1e6);
    end

    #=
    # plot the results
    figure(1)
    clf()
    plot(t ./ 3600, Agpp_t, "g-");
    plot(t ./ 3600, Asum_t, "k-");

    figure(2)
    clf()
    plot(t ./ 3600, Esum_t, "k-");
    =#

    return nothing
end
