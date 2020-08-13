###############################################################################
#
# Calculate SIF fluxes
#
###############################################################################
"""
    sif_fluxes!(leaf_array::Array{LeafBios{FT},1}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, can::Canopy4RT{FT}, soil_opt::SoilOpticals{FT}, wl_set::WaveLengths{FT}) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport for SIF radiation (calls `compute_diffusive_S` internally).
Layer reflectance and transmission is computed from SW optical properties, layer sources from absorbed light and SIF efficiencies. Boundary conditions are zero SIF incoming from atmosphere or soil.
- `leaf_array` An array of [`LeafBios`](@ref) type struct (i.e. leaf optical properties can change with canopy height)
- `can_opt` A [`CanopyOpticals`](@ref) struct for providing optical layer properties
- `can_rad` A [`CanopyRads`](@ref) struct
- `can` A [`Canopy4RT`](@ref) type struct for providing LAI and nLayer and clumping
- `soil_opt` A [`SoilOpticals`](@ref) type struct for soil optical properties
- `wl_set` An [`WaveLengths`](@ref) type struct
"""
function sif_fluxes!(
            leaf_array::Array{LeafBios{FT},1},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            can::Canopy4RT{FT},
            soil_opt::SoilOpticals{FT},
            wl_set::WaveLengths{FT},
            rt_con::RTContainer{FT}
) where {FT<:AbstractFloat}
    # 1. unpack variables from structures
    @unpack LAI, lidf, nLayer, Ω = can;
    @unpack a, absfo, absfs, absfsfo, cosΘ_l, cos2Θ_l, fo, fs, fsfo, Po, Ps,
            Pso, sigb, vb, vf = can_opt;
    @unpack E_down, E_up,ϕ_shade,ϕ_sun = can_rad;
    @unpack soil_sif_albedo = rt_con;
    @unpack dWL, iWLE, iWLF, nWLF = wl_set;

    # 2. calculate some useful parameters
    iLAI = Ω * LAI / nLayer;
    rt_con.τ_dd_sif .= 1 .- view(a, iWLF, :) .* iLAI;
    rt_con.ρ_dd_sif .= view(sigb, iWLF, :) .* iLAI;
    @unpack τ_dd_sif, ρ_dd_sif = rt_con;

    # 3. Compute mid layer Ps,Po,Pso
    #    Qso = (Pso[1:nLayer] + Pso[2:nLayer+1]) / 2;
    #    Qs  = ( Ps[1:nLayer] +  Ps[2:nLayer+1]) / 2;
    #    Qo  = ( Po[1:nLayer] +  Po[2:nLayer+1]) / 2;
    Qso = view(Pso, 1:nLayer);
    Qs  = view(Ps , 1:nLayer);
    Qo  = view(Po , 1:nLayer);

    # 4.  reflectance, transmittance factors in a thin layer the following are
    #     vectors with length [nl,nWL]
    rt_con.sun_dwl_iWlE .= view(can_opt.Es_, iWLE, 1) .* rt_con.dλ_iWlE;
    @inbounds for i=1:nLayer
        if length(leaf_array)>1
            Mb = leaf_array[i].Mb;
            Mf = leaf_array[i].Mf;
        else
            Mb = leaf_array[1].Mb;
            Mf = leaf_array[1].Mf;
        end
        rt_con.M⁺ .= (Mb .+ Mf) ./ 2;
        rt_con.M⁻ .= (Mb .- Mf) ./ 2;

        # Need to normalize incoming radiation bin
        # change from mSCOPE to enable different WL grids!
        mul!(rt_con.M⁻_sun, rt_con.M⁻, rt_con.sun_dwl_iWlE);
        mul!(rt_con.M⁺_sun, rt_con.M⁺, rt_con.sun_dwl_iWlE);

        rt_con.tmp_dwl_iWlE .= view(E_down, iWLE, i) .* rt_con.dλ_iWlE;
        mul!(rt_con.M⁺⁻, rt_con.M⁺, rt_con.tmp_dwl_iWlE);
        rt_con.tmp_dwl_iWlE .= view(E_up, iWLE, i+1) .* rt_con.dλ_iWlE;
        mul!(rt_con.M⁺⁺, rt_con.M⁺, rt_con.tmp_dwl_iWlE);
        rt_con.tmp_dwl_iWlE .= view(E_up, iWLE, i+1) .* rt_con.dλ_iWlE;
        mul!(rt_con.M⁻⁺, rt_con.M⁻, rt_con.tmp_dwl_iWlE);
        rt_con.tmp_dwl_iWlE .= view(E_down, iWLE, i) .* rt_con.dλ_iWlE;
        mul!(rt_con.M⁻⁻, rt_con.M⁻, rt_con.tmp_dwl_iWlE);

        # Here comes the tedious part:
        # TODO move them to a seprate function
        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* cosΘ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        sunCos = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= ϕ_shade[i] .* cosΘ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        shadeCos = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* cos2Θ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        sunCos2 = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= ϕ_shade[i] .* cos2Θ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        shadeCos2 = mean(rt_con.ϕ_cosΘ_lidf);
        mul!(rt_con.ϕ_cosΘ_lidf, view(ϕ_sun, :, :, i)', lidf);
        sunLidf = mean(rt_con.ϕ_cosΘ_lidf);
        shadeLidf = mean(lidf) * ϕ_shade[i];

        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfsfo;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_absfsfo = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fsfo;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_fsfo = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.wfEs .= _mean_absfsfo .* rt_con.M⁺_sun .+
                       _mean_fsfo    .* rt_con.M⁻_sun;

        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfs;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_absfs = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fs .* cosΘ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_fs = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.sfEs .= _mean_absfs .* rt_con.M⁺_sun .-
                       _mean_fs    .* rt_con.M⁻_sun;
        rt_con.sbEs .= _mean_absfs .* rt_con.M⁺_sun .+
                       _mean_fs    .* rt_con.M⁻_sun;

        rt_con.ϕ_cosΘ .= ϕ_shade[i] .* absfo;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_absfo = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= ϕ_shade[i] .* fo .* cosΘ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_fo = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.vfEplu_shade .= _mean_absfo .* rt_con.M⁺⁺ .-
                               _mean_fo    .* rt_con.M⁻⁺;
        rt_con.vbEmin_shade .= _mean_absfo .* rt_con.M⁺⁻ .+
                               _mean_fo    .* rt_con.M⁻⁻;

        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfo;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_absfo = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fo .* cosΘ_l;
        mul!(rt_con.ϕ_cosΘ_lidf, (rt_con.ϕ_cosΘ)', lidf);
        _mean_fo = mean(rt_con.ϕ_cosΘ_lidf);
        rt_con.vfEplu_sun .= _mean_absfo .* rt_con.M⁺⁺ .- _mean_fo .* rt_con.M⁻⁺;
        rt_con.vbEmin_sun .= _mean_absfo .* rt_con.M⁺⁻ .+ _mean_fo .* rt_con.M⁻⁻;

        rt_con.sigfEmin_shade .= shadeLidf .* rt_con.M⁺⁻ .-
                                 shadeCos2 .* rt_con.M⁻⁻;
        rt_con.sigbEmin_shade .= shadeLidf .* rt_con.M⁺⁻ .+
                                 shadeCos2 .* rt_con.M⁻⁻;
        rt_con.sigfEmin_sun   .= sunLidf   .* rt_con.M⁺⁻ .-
                                 sunCos2   .* rt_con.M⁻⁻;
        rt_con.sigbEmin_sun   .= sunLidf   .* rt_con.M⁺⁻ .+
                                 sunCos2   .* rt_con.M⁻⁻;
        rt_con.sigfEplu_shade .= shadeLidf .* rt_con.M⁺⁺ .-
                                 shadeCos2 .* rt_con.M⁻⁺;
        rt_con.sigbEplu_shade .= shadeLidf .* rt_con.M⁺⁺ .+
                                 shadeCos2 .* rt_con.M⁻⁺;
        rt_con.sigfEplu_sun   .= sunLidf   .* rt_con.M⁺⁺ .-
                                 sunCos2   .* rt_con.M⁻⁺;
        rt_con.sigbEplu_sun   .= sunLidf   .* rt_con.M⁺⁺ .+
                                 sunCos2   .* rt_con.M⁻⁺;

        # Fluxes:
        #    sunlit for each layer
        #    shade leaf for each layer
        #    Eq. 29a for sunlit leaf
        #    Eq. 29b for sunlit leaf
        #    Eq. 29a for shade leaf
        #    Eq. 29b for shade leaf
        rt_con.piLs[ :,i] .= rt_con.wfEs           .+
                             rt_con.vfEplu_sun     .+
                             rt_con.vbEmin_sun;
        rt_con.piLd[ :,i] .= rt_con.vbEmin_shade   .+
                             rt_con.vfEplu_shade;
        rt_con.Fsmin[:,i] .= rt_con.sfEs           .+
                             rt_con.sigfEmin_sun   .+
                             rt_con.sigbEplu_sun;
        rt_con.Fsplu[:,i] .= rt_con.sbEs           .+
                             rt_con.sigbEmin_sun   .+
                             rt_con.sigfEplu_sun;
        rt_con.Fdmin[:,i] .= rt_con.sigfEmin_shade .+
                             rt_con.sigbEplu_shade;
        rt_con.Fdplu[:,i] .= rt_con.sigbEmin_shade .+
                             rt_con.sigfEplu_shade;





        # Total weighted fluxes
        _qs_iLAI   = Qs[i] * iLAI;
        _1_qs_iLAI = (1 - Qs[i]) * iLAI;
        rt_con.S⁻[:,i]   .= _qs_iLAI   .* view(rt_con.Fsmin, :, i) .+
                            _1_qs_iLAI .* view(rt_con.Fdmin, :, i);
        rt_con.S⁺[:,i]   .= _qs_iLAI   .* view(rt_con.Fsplu, :, i) .+
                            _1_qs_iLAI .* view(rt_con.Fdplu, :, i);
        rt_con.Femo[:,i] .= _qs_iLAI   .* view(rt_con.piLs , :, i) .+
                            _1_qs_iLAI .* view(rt_con.piLd , :, i);
    end

    # 5. Compute diffusive fluxes within canopy
    #    Use Zero SIF fluxes as top and bottom boundary:
    # TODO pre-allocate these!
    diffusive_S!(rt_con.F⁻,
                 rt_con.F⁺,
                 rt_con.net_diffuse,
                 τ_dd_sif,
                 ρ_dd_sif,
                 rt_con.S⁻,
                 rt_con.S⁺,
                 rt_con.zeroB,
                 rt_con.zeroB,
                 soil_sif_albedo);



    # 6. Save in output structures!
    #    direct Sunlit leaves
    #    direct shaded leaves
    _iLAI_pi = iLAI / FT(pi);

    # why this step is so slow?
    rt_con.tmp_1d_nLayer .= Qso .* _iLAI_pi;
    mul!(can_rad.SIF_obs_sunlit, rt_con.piLs, rt_con.tmp_1d_nLayer);
    rt_con.tmp_1d_nLayer .= (Qo .- Qso) .* _iLAI_pi;
    mul!(can_rad.SIF_obs_shaded, rt_con.piLd, rt_con.tmp_1d_nLayer);

    # 7. SIF from scattered internally and soil contribution
    rt_con.tmp_2d_nWlF_nLayer .= view(vb, iWLF, :) .* view(rt_con.F⁻, :, 1:nLayer) .+
                                 view(vf, iWLF, :) .* view(rt_con.F⁺, :, 1:nLayer);
    mul!(can_rad.SIF_obs_scattered, rt_con.tmp_2d_nWlF_nLayer, Qo);
    can_rad.SIF_obs_scattered .= can_rad.SIF_obs_scattered .* _iLAI_pi;
    can_rad.SIF_obs_soil .= ( soil_sif_albedo .* view(rt_con.F⁻, :, nLayer+1) ) .* Po[end] ./ FT(pi);

    can_rad.SIF_hemi .= view(rt_con.F⁺, :, 1);
    can_rad.SIF_obs  .= can_rad.SIF_obs_sunlit    .+
                        can_rad.SIF_obs_shaded    .+
                        can_rad.SIF_obs_scattered .+
                        can_rad.SIF_obs_soil;
    # can_rad.SIF_sum[:]  = sum(rt_con.S⁻ + rt_con.S⁺, dims=2);
    @inbounds for j in eachindex(can_rad.SIF_sum)
        can_rad.SIF_sum[j] = sum( view(rt_con.S⁻, j, :) ) +
                             sum( view(rt_con.S⁺, j, :) );
    end

    return nothing
end
