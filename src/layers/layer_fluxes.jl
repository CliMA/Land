###############################################################################
#
# Update g_sw using different stomatal models
# This part should be move to StomataModels in its next release
#
###############################################################################
"""
    update_gsw!(clayer::CanopyLayer{FT},
                sm::EmpiricalStomatalModel{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                envir::AirLayer{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}
    update_gsw!(clayer::CanopyLayer{FT},
                sm::OSMWang{FT},
                photo_set::AbstractPhotoModelParaSet{FT},
                envir::AirLayer{FT},
                Δt::FT
    ) where {FT<:AbstractFloat}

Update g_sw prognostically, given
- `clayer` A canoy layer --- `CanopyLayer` type struct
- `sm` `EmpiricalStomatalModel` or `OSMWang` type stomatal model
- `photo_set` AbstractPhotoModelParaSet type photosynthesis model
- `envir` `AirLayer` type environmental conditions
- `Δt` Time interval for prognostic stomatal conductance
"""
function update_gsw!(
            clayer::CanopyLayer{FT},
            sm::EmpiricalStomatalModel{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    # unpack values
    @unpack g_sw, n_leaf, ps = clayer;

    # calculate steady state values
    #
    #
    #
    #
    # use a container here in CanopyLayer in the release version
    #
    #
    #
    #
    gsw_ss = max.(sm.g0, empirical_gsw_from_model(sm, clayer, envir, FT(1)));
    # assume τ = 10 minutes
    for iLF in 1:n_leaf
        g_sw[iLF] += (gsw_ss[iLF] - g_sw[iLF]) / 600 * Δt;
    end

    return nothing
end




function update_gsw!(
            clayer::CanopyLayer{FT},
            sm::OSMWang{FT},
            photo_set::AbstractPhotoModelParaSet{FT},
            envir::AirLayer{FT},
            Δt::FT
) where {FT<:AbstractFloat}
    # unpack values
    @unpack An, ec, g_bc, g_bw, g_lw, g_m, g_sw, n_leaf, p_sat, ps = clayer;
    @unpack p_atm, p_H₂O = envir;

    # update g_sw
    for iLF in 1:n_leaf
        _gsw = g_sw[iLF] .+ FT(0.001);
        _glw = 1 / ( 1/_gsw + 1/g_bw[iLF] );
        _gsc = _gsw / FT(1.6);
        _glc = 1 / ( 1/_gsc +  1/g_m[iLF] +  1/g_bc[iLF] );
        leaf_photo_from_glc!(photo_set, ps, envir, _glc);
        ∂A   = ps.An - An[iLF];
        e0   = g_lw[iLF] * (p_sat - p_H₂O) / p_atm;
        e1   = _glw * (p_sat - p_H₂O) / p_atm;
        ∂E   = e1 - e0;
        ∂A∂E = ∂A / ∂E;
        ∂Θ∂E = FT(max(0.1, An[iLF]) / max(ec - e0, 1e-7));
        #
        #
        #
        #
        # remember to add a τ parameter in CanopyLayer
        # assume τ is 1e-6 here
        #
        #
        #
        #
        g_sw[iLF] += (∂A∂E - ∂Θ∂E) * FT(1e-6) * Δt;
    end

    return nothing
end








###############################################################################
#
# Update carbon and water fluxes
#
###############################################################################
"""
    layer_fluxes!(node::SPACMono{FT}) where {FT<:AbstractFloat}

Run carbon, water, and energy fluxes for all canopy layers, given
- `node` [`SPACMono`](@ref) type struct
"""
function layer_fluxes!(node::SPACMono{FT}) where {FT<:AbstractFloat}
    # 0.1 unpack data
    @unpack angles, can_opt, can_rad, canopy_rt, envirs, f_SL, ga, in_rad,
            leaves_rt, n_canopy, photo_set, plant_hs, plant_ps, rt_con,
            soil_opt, stomata_model, wl_set = node;
    @unpack nAzi, nIncl = canopy_rt;
    nSL = nAzi * nIncl;

    canopy_geometry!(canopy_rt, angles, can_opt, rt_con);
    canopy_matrices!(leaves_rt, can_opt);
    short_wave!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, rt_con);
    canopy_fluxes!(canopy_rt, can_opt, can_rad, in_rad, soil_opt, leaves_rt,
                   wl_set, rt_con);

    #
    f_H₂O = 0;
    f_GPP = 0;
    f_NPP = 0;
    for i_can in 1:n_canopy
        iEN = envirs[i_can];
        iHS = plant_hs.leaves[i_can];
        iPS = plant_ps[i_can];
        iRT = n_canopy + 1 - i_can;

        iPS.T = can_rad.T_sun[i_can];
        update_leaf_TP!(photo_set, iPS, iHS, iEN);
        temperature_effects!(iHS, iPS.T);

        # calculate the fraction of sunlit and shaded leaves
        f_view = (can_opt.Ps[iRT] + can_opt.Ps[iRT+1]) / 2;
        for iLF in 1:nSL
            iPS.APAR[iLF] = can_rad.absPAR_sunCab[(iRT-1)*nSL+iLF] * FT(1e6);
            iPS.LAIx[iLF] = f_view * f_SL[iLF];
        end
        iPS.APAR[end] = can_rad.absPAR_shadeCab[iRT] * FT(1e6);
        iPS.LAIx[end] = 1 - f_view;

        # iterate for 15 times to find steady state solution
        for iter in 1:15
            # calculate the photosynthetic rates
            update_leaf_from_gsw!(photo_set, iPS, iEN);
            update_gsw!(iPS, stomata_model, photo_set, iEN, FT(120));
            leaf_gsw_control!(photo_set, iPS, iEN);
        end

        # update the flow rates
        for iLF in 1:(nSL+1)
            f_GPP += iPS.Ag[iLF] * iPS.LAIx[iLF] * iPS.LA;
            f_NPP += iPS.An[iLF] * iPS.LAIx[iLF] * iPS.LA;
            f_H₂O += iPS.g_lw[iLF] * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm *
                     iPS.LAIx[iLF] * iPS.LA;
        end
    end

    # update the flows in SPACMono
    node.f_gpp = f_GPP / ga;
    node.f_npp = f_NPP / ga;
    node.f_H₂O = f_H₂O / ga;

    return nothing
end
