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

function test_diurnal(
            node::SPACMono{FT},
            Δt::FT=FT(10)
) where {FT<:AbstractFloat}
    # 0.1 create variables required
    @unpack n_canopy, wl_set = node;
    canopy_rt = node.canopy_rt;
    fraction_sl::Array{FT,1} = repeat(canopy_rt.lidf, outer=[ length(canopy_rt.lazitab) ]) / length(canopy_rt.lazitab);
    n_sl = length(canopy_rt.lidf) * length(canopy_rt.lazitab);




    # change node.in_rad to alter direct and diffuse fraction
    # sum of direct  energy in W:
    # e_dir = sum( (wl_set.swl[2:end] .- wl_set.swl[1:end-1]) .* node.in_rad.E_direct  ) / 1000
    # sum of diffuse energy in W:
    # e_dif = sum( (wl_set.swl[2:end] .- wl_set.swl[1:end-1]) .* node.in_rad.E_diffuse ) / 1000
    # target e_dir is tar_dir, do it this way:
    # node.in_rad.E_direct  .*= tar_dir / e_dir
    # target e_dif is tar_dif, do it this way:
    # node.in_rad.E_diffuse .*= tar_dif / e_dif



    # 0.2 A diurnal cycle for radiation and Tair
    Tmean  = FT(295.15)
    DeltaT = 3
    omega  = 2 * π / (24*3600)
    t      = range(0, stop=24*3600*2, step=Δt)
    phi_t  = omega .* t .- π .* ones(FT, size(t)) ./ 2
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
        leaf_temperature_dependence!(node.photo_set, node.plant_ps[i_can].ps, node.envirs[i_can], Tlef_t[1]);
        node.envirs[i_can].p_H₂O = node.envirs[i_can].p_sat / 2;
        # latent heat flux as well
        node.plant_ps[i_can].LV = latent_heat_vapor(Tlef_t[i_can]) * 1000 / 18;
        # vc dependecy as well
        temperature_effects!(node.plant_hs.leaves[i_can], Tlef_t[i_can]);
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
        zenith = zenith_angle(node.latitude, 180 + t[i_tim]/86400 );
        zenith = min(88, zenith);
        node.angles.tts = zenith;
        canopy_geometry!(canopy_rt, node.angles, node.can_opt, node.rt_con);
        canopy_matrices!(node.leaves_rt, node.can_opt);
        short_wave!(canopy_rt, node.can_opt, node.can_rad, node.in_rad, node.soil_opt, node.rt_con);
        canopy_fluxes!(canopy_rt, node.can_opt, node.can_rad, node.in_rad, node.soil_opt, node.leaves_rt, wl_set, node.rt_con);
        thermal_fluxes!(node.leaves_rt, node.can_opt, node.can_rad, canopy_rt, node.soil_opt, [FT(400.0)], wl_set);
        SIF_fluxes!(node.leaves_rt, node.can_opt, node.can_rad, canopy_rt, node.soil_opt, wl_set, node.rt_con, node.rt_dim);

        # 2. Update transpiration rate
        a_gpp = FT(0);
        a_sum = FT(0);
        e_sum = FT(0);
        for i_can in 1:n_canopy
            rt_layer = node.n_canopy + 1 - i_can;

            # set the diffuse par to all the leaves
            #node.plant_ps[i_can].APAR           .= PARr_t[i_tim]  * node.can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
            #node.plant_ps[i_can].APAR[1:end-1] .+= PARr_t[i_tim] .* reshape(node.can_rad.absPAR_sunCab[:,:,rt_layer],(:,1))[:,1] .* FT(1e6);

            # calculate the fraction of sunlit and shaded leaves
            f_view = (node.can_opt.Ps[rt_layer] + node.can_opt.Ps[rt_layer+1]) / 2;
            #la_new = [f_view .* fraction_sl; 1-f_view];
            #node.plant_ps[i_can].LAIx .= la_new;
            ##=
            for i_leaf in 1:n_sl
                node.plant_ps[i_can].APAR[i_leaf] = PARr_t[i_tim] * node.can_rad.absPAR_sunCab[(rt_layer-1)*n_sl + i_leaf] * FT(1e6);
                node.plant_ps[i_can].LAIx[i_leaf] = f_view * fraction_sl[i_leaf];
            end
            node.plant_ps[i_can].APAR[end] = PARr_t[i_tim] * node.can_rad.absPAR_shadeCab[rt_layer] * FT(1e6);
            node.plant_ps[i_can].LAIx[end] = 1 - f_view;
            # =#

            # update leaf temperature
            # update the parameters from stomatal model
            node.plant_ps[i_can].T = Tlef_t[i_tim];
            update_leaf_TP!(node.photo_set, node.plant_ps[i_can], node.plant_hs.leaves[i_can], node.envirs[i_can]);
            # uncomment this for Sperry, Eller, and Gentine models
            #update_leaf_AK!(node.photo_set, node.plant_ps[i_can], node.plant_hs.leaves[i_can], node.envirs[i_can]);

            #canopyi.t_list[1:end-1] .= reshape(node.can_rad.T_sun3D[:,:,rt_layer],(:,1))[:,1]
            #canopyi.t_list[end]      = node.can_rad.T_shade[rt_layer]
            # Assume RH = 50%
            node.envirs[i_can].t_air = Tair_t[i_tim];
            node.envirs[i_can].p_sat = saturation_vapor_pressure(Tair_t[i_tim]);
            node.envirs[i_can].p_H₂O = node.envirs[i_can].p_sat / 2;

            # update the flow rates
            g_i_can = FT(0);
            a_i_can = FT(0);
            e_i_can = FT(0);
            for i_leaf in 1:(n_sl+1)
                g_i_can += node.plant_ps[i_can].Ag[i_leaf] *
                           node.plant_ps[i_can].LAIx[i_leaf] *
                           node.plant_ps[i_can].LA;
                a_i_can += node.plant_ps[i_can].An[i_leaf] *
                           node.plant_ps[i_can].LAIx[i_leaf] *
                           node.plant_ps[i_can].LA;
                e_i_can += node.plant_ps[i_can].g_lw[i_leaf] *
                           (node.plant_ps[i_can].p_sat - node.envirs[i_can].p_H₂O) /
                           node.envirs[i_can].p_atm *
                           node.plant_ps[i_can].LAIx[i_leaf] *
                           node.plant_ps[i_can].LA;
            end
            #g_i_can = sum( node.plant_ps[i_can].Ag .* node.plant_ps[i_can].LAIx ) * node.plant_ps[i_can].LA;
            #a_i_can = sum( node.plant_ps[i_can].An .* node.plant_ps[i_can].LAIx ) * node.plant_ps[i_can].LA;
            #e_i_can = sum( node.plant_ps[i_can].g_lw .* (node.plant_ps[i_can].p_sat - node.envirs[i_can].p_H₂O) ./ node.envirs[i_can].p_atm .* node.plant_ps[i_can].LAIx ) * node.plant_ps[i_can].LA;

            # calculate the photosynthetic rates
            gas_exchange!(node.photo_set, node.plant_ps[i_can], node.envirs[i_can], GswDrive());
            gsw_control!(node.photo_set, node.plant_ps[i_can], node.envirs[i_can]);

            # use the ball-berry model here for now as the ∂A/∂E and ∂A/∂Θ functions are not yet ready
            gsw_ss = stomatal_conductance(node.stomata_model, node.plant_ps[i_can], node.envirs[i_can], FT(1));

            # assume τ = 10 minutes
            for i_leaf in 1:(n_sl+1)
                node.plant_ps[i_can].g_sw[i_leaf] += (gsw_ss[i_leaf] - node.plant_ps[i_can].g_sw[i_leaf]) / 600 * Δt;
            end
            #node.plant_ps[i_can].g_sw .+= (gsw_ss .- node.plant_ps[i_can].g_sw) / 600 * Δt;

            a_gpp += g_i_can;
            a_sum += a_i_can;
            e_sum += e_i_can;
        end

        Agpp_t[i_tim] = a_gpp / node.ga;
        Asum_t[i_tim] = a_sum / node.ga;
        Esum_t[i_tim] = e_sum / node.ga;

        # update pressure profile every 10 minutes
        # TO BE DONE within PlantHydraulics

        _par = PARr_t[i_tim]  * node.can_rad.absPAR_shadeCab[1] * FT(1e6);
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
