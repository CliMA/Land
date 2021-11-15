###############################################################################
#
# Calculate SIF fluxes
#
###############################################################################
"""
    SIF_fluxes!(leaves::Array{LeafBios{FT},1},
                can_opt::CanopyOpticals{FT},
                can_rad::CanopyRads{FT},
                can::Canopy4RT{FT},
                soil::SoilOpticals{FT},
                wls::WaveLengths{FT},
                rt_con::RTCache{FT},
                rt_dim::RTDimensions;
                photon::Bool = true
    ) where {FT<:AbstractFloat}

Computes 2-stream diffusive radiation transport for SIF radiation (calls
    [`diffusive_S!`] internally). Layer reflectance and transmission is
    computed from SW optical properties, layer sources from absorbed light and
    SIF efficiencies. Boundary conditions are zero SIF incoming from atmosphere
    or soil.
- `leaves` Array of [`LeafBios`](@ref) type struct
- `can_opt` [`CanopyOpticals`](@ref) type struct
- `can_rad` [`CanopyRads`](@ref) type struct
- `can` [`Canopy4RT`](@ref) type struct
- `soil` [`SoilOpticals`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
- `rt_dim` [`RTDimensions`](@ref) type struct
- `photon` If true, use photon unit in the matrix conversion

"""
function SIF_fluxes!(
            leaves::Array{LeafBios{FT},1},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            can::Canopy4RT{FT},
            soil::SoilOpticals{FT},
            wls::WaveLengths{FT},
            rt_con::RTCache{FT},
            rt_dim::RTDimensions;
            photon::Bool = true
) where {FT<:AbstractFloat}
    # unpack variables from structures
    @unpack LAI, lidf, nLayer, Ω = can;
    @unpack a, absfo, absfs, absfsfo, cosΘ_l, cos2Θ_l, fo, fs, fsfo, Po, Ps, Pso, sigb, vb, vf = can_opt;
    @unpack E_down, E_up, ϕ_shade, ϕ_sun = can_rad;
    @unpack ρ_SW_SIF = soil;
    @unpack dWL_iWLE, iWLE, iWLF, WLE, WLF = wls;
    sf_con = rt_con.sf_con;

    # 1. define some useful parameters
    iLAI = LAI * Ω / nLayer;

    # 2. calculate some useful parameters
    sf_con.τ_dd .= 1 .- view(a, iWLF, :) .* iLAI;
    sf_con.ρ_dd .= view(sigb, iWLF, :) .* iLAI;

    # 3. Compute mid layer Ps,Po,Pso
    #    Qso = (Pso[1:nLayer] + Pso[2:nLayer+1]) / 2;
    #    Qs  = ( Ps[1:nLayer] +  Ps[2:nLayer+1]) / 2;
    #    Qo  = ( Po[1:nLayer] +  Po[2:nLayer+1]) / 2;
    Qso = view(Pso, 1:nLayer);
    Qs  = view(Ps , 1:nLayer);
    Qo  = view(Po , 1:nLayer);

    # 4.  reflectance, transmittance factors in a thin layer the following are
    #     vectors with length [nl,nWL]
    sf_con.sun_dwl_iWlE .= view(can_opt.Es_, iWLE, 1) .* dWL_iWLE;
    if photon sf_con.sun_dwl_iWlE .*= WLE .* _FAC(FT) end;
    @inbounds for i=1:nLayer
        if length(leaves)>1
            Mb = leaves[i].Mb;
            Mf = leaves[i].Mf;
        else
            Mb = leaves[1].Mb;
            Mf = leaves[1].Mf;
        end
        sf_con.M⁺ .= (Mb .+ Mf) ./ 2;
        sf_con.M⁻ .= (Mb .- Mf) ./ 2;

        # Need to normalize incoming radiation bin
        # change from mSCOPE to enable different WL grids!
        mul!(sf_con.M⁻_sun, sf_con.M⁻, sf_con.sun_dwl_iWlE);
        mul!(sf_con.M⁺_sun, sf_con.M⁺, sf_con.sun_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_down, iWLE, i) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁺⁻, sf_con.M⁺, sf_con.tmp_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_up, iWLE, i+1) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁺⁺, sf_con.M⁺, sf_con.tmp_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_up, iWLE, i+1) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁻⁺, sf_con.M⁻, sf_con.tmp_dwl_iWlE);

        sf_con.tmp_dwl_iWlE .= view(E_down, iWLE, i) .* dWL_iWLE;
        if photon sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT) end;
        mul!(sf_con.M⁻⁻, sf_con.M⁻, sf_con.tmp_dwl_iWlE);

        # Here comes the tedious part:
        # TODO move them to a seprate function
        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        sunCos = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        shadeCos = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* cos2Θ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        sunCos2 = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* cos2Θ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        shadeCos2 = mean(sf_con.ϕ_cosΘ_lidf);

        mul!(sf_con.ϕ_cosΘ_lidf, view(ϕ_sun, :, :, i)', lidf);
        sunLidf = mean(sf_con.ϕ_cosΘ_lidf);
        shadeLidf = mean(lidf) * ϕ_shade[i];

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfsfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfsfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fsfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fsfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.wfEs .= _mean_absfsfo .* sf_con.M⁺_sun .+
                       _mean_fsfo    .* sf_con.M⁻_sun;

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfs;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfs = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fs .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fs = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.sfEs .= _mean_absfs .* sf_con.M⁺_sun .- _mean_fs .* sf_con.M⁻_sun;
        sf_con.sbEs .= _mean_absfs .* sf_con.M⁺_sun .+ _mean_fs .* sf_con.M⁻_sun;

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* absfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= ϕ_shade[i] .* fo .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.vfEplu_shade .= _mean_absfo .* sf_con.M⁺⁺ .- _mean_fo .* sf_con.M⁻⁺;
        sf_con.vbEmin_shade .= _mean_absfo .* sf_con.M⁺⁻ .+ _mean_fo .* sf_con.M⁻⁻;

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* absfo;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_absfo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.ϕ_cosΘ .= view(ϕ_sun, :, :, i) .* fo .* cosΘ_l;
        mul!(sf_con.ϕ_cosΘ_lidf, adjoint(sf_con.ϕ_cosΘ), lidf);
        _mean_fo = mean(sf_con.ϕ_cosΘ_lidf);

        sf_con.vfEplu_sun .= _mean_absfo .* sf_con.M⁺⁺ .- _mean_fo .* sf_con.M⁻⁺;
        sf_con.vbEmin_sun .= _mean_absfo .* sf_con.M⁺⁻ .+ _mean_fo .* sf_con.M⁻⁻;

        sf_con.sigfEmin_shade .= shadeLidf .* sf_con.M⁺⁻ .- shadeCos2 .* sf_con.M⁻⁻;
        sf_con.sigbEmin_shade .= shadeLidf .* sf_con.M⁺⁻ .+ shadeCos2 .* sf_con.M⁻⁻;
        sf_con.sigfEmin_sun   .= sunLidf   .* sf_con.M⁺⁻ .- sunCos2   .* sf_con.M⁻⁻;
        sf_con.sigbEmin_sun   .= sunLidf   .* sf_con.M⁺⁻ .+ sunCos2   .* sf_con.M⁻⁻;
        sf_con.sigfEplu_shade .= shadeLidf .* sf_con.M⁺⁺ .- shadeCos2 .* sf_con.M⁻⁺;
        sf_con.sigbEplu_shade .= shadeLidf .* sf_con.M⁺⁺ .+ shadeCos2 .* sf_con.M⁻⁺;
        sf_con.sigfEplu_sun   .= sunLidf   .* sf_con.M⁺⁺ .- sunCos2   .* sf_con.M⁻⁺;
        sf_con.sigbEplu_sun   .= sunLidf   .* sf_con.M⁺⁺ .+ sunCos2   .* sf_con.M⁻⁺;

        # Fluxes:
        #    sunlit for each layer
        #    shade leaf for each layer
        #    Eq. 29a for sunlit leaf
        #    Eq. 29b for sunlit leaf
        #    Eq. 29a for shade leaf
        #    Eq. 29b for shade leaf
        sf_con.piLs[:,i]  .= sf_con.wfEs .+ sf_con.vfEplu_sun .+ sf_con.vbEmin_sun;
        sf_con.piLd[:,i]  .= sf_con.vbEmin_shade .+ sf_con.vfEplu_shade;
        sf_con.Fsmin[:,i] .= sf_con.sfEs .+ sf_con.sigfEmin_sun .+ sf_con.sigbEplu_sun;
        sf_con.Fsplu[:,i] .= sf_con.sbEs .+ sf_con.sigbEmin_sun .+ sf_con.sigfEplu_sun;
        sf_con.Fdmin[:,i] .= sf_con.sigfEmin_shade .+ sf_con.sigbEplu_shade;
        sf_con.Fdplu[:,i] .= sf_con.sigbEmin_shade .+ sf_con.sigfEplu_shade;

        # Total weighted fluxes
        _qs_iLAI   = Qs[i] * iLAI;
        _1_qs_iLAI = (1 - Qs[i]) * iLAI;
        sf_con.S⁻[:,i]   .= _qs_iLAI .* view(sf_con.Fsmin, :, i) .+ _1_qs_iLAI .* view(sf_con.Fdmin, :, i);
        sf_con.S⁺[:,i]   .= _qs_iLAI .* view(sf_con.Fsplu, :, i) .+ _1_qs_iLAI .* view(sf_con.Fdplu, :, i);
        sf_con.Femo[:,i] .= _qs_iLAI .* view(sf_con.piLs , :, i) .+ _1_qs_iLAI .* view(sf_con.piLd , :, i);
    end

    # 5. Compute diffusive fluxes within canopy
    #    Use Zero SIF fluxes as top and bottom boundary:
    # TODO pre-allocate these!
    diffusive_S!(sf_con, soil, rt_dim);

    # 6. Save in output structures!
    #    direct Sunlit leaves
    #    direct shaded leaves
    _iLAI_pi = iLAI / FT(pi);

    # why this step is so slow?
    sf_con.tmp_1d_nLayer .= Qso .* _iLAI_pi;
    mul!(can_rad.SIF_obs_sunlit, sf_con.piLs, sf_con.tmp_1d_nLayer);
    sf_con.tmp_1d_nLayer .= (Qo .- Qso) .* _iLAI_pi;
    mul!(can_rad.SIF_obs_shaded, sf_con.piLd, sf_con.tmp_1d_nLayer);

    # 7. SIF from scattered internally and soil contribution
    sf_con.tmp_2d_nWlF_nLayer .= view(vb, iWLF, :) .* view(sf_con.F⁻, :, 1:nLayer) .+ view(vf, iWLF, :) .* view(sf_con.F⁺, :, 1:nLayer);
    mul!(can_rad.SIF_obs_scattered, sf_con.tmp_2d_nWlF_nLayer, Qo);
    can_rad.SIF_obs_scattered .= can_rad.SIF_obs_scattered .* _iLAI_pi;
    can_rad.SIF_obs_soil .= ( ρ_SW_SIF .* view(sf_con.F⁻, :, nLayer+1) ) .* Po[end] ./ FT(pi);

    can_rad.SIF_hemi .= view(sf_con.F⁺, :, 1);
    can_rad.SIF_obs  .= can_rad.SIF_obs_sunlit .+ can_rad.SIF_obs_shaded .+ can_rad.SIF_obs_scattered .+ can_rad.SIF_obs_soil;
    # can_rad.SIF_sum[:]  = sum(sf_con.S⁻ + sf_con.S⁺, dims=2);
    @inbounds for j in eachindex(can_rad.SIF_sum)
        can_rad.SIF_sum[j] = sum( view(sf_con.S⁻, j, :) ) + sum( view(sf_con.S⁺, j, :) );
    end

    if photon
        can_rad.SIF_obs ./= WLF .* _FAC(FT);
    end;

    return nothing
end




"""
    SIF_fluxes!(leaf::LeafBios{FT},
                in_rad::IncomingRadiation{FT},
                wls::WaveLengths{FT},
                rt_con::RTCache{FT},
                fqe::FT = FT(0.01);
                photon::Bool = true
    ) where {FT<:AbstractFloat}

Leaf level SIF, given
- `leaf` [`LeafBios`](@ref) type struct
- `in_rad` [`IncomingRadiation`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache
- `fqe` Fluorescence quantum yield (default at 1%)
- `photon` If true, use photon unit in the matrix conversion

Note that `in_rad` assumes direct light with zenith angle of 0, and a zenith
    angle correction needs to be made before passing it to this function. The
    up- and down-ward SIF are stored in `sf_con` as `M⁻_sun` and `M⁺_sun`.
"""
function SIF_fluxes!(
            leaf::LeafBios{FT},
            in_rad::IncomingRadiation{FT},
            wls::WaveLengths{FT},
            rt_con::RTCache{FT},
            fqe::FT = FT(0.01);
            photon::Bool = true
) where {FT<:AbstractFloat}
    # unpack the values
    @unpack Mb, Mf = leaf;
    @unpack dWL_iWLE, iWLE, WLE, WLF = wls;
    sf_con = rt_con.sf_con;
    sf_con.tmp_dwl_iWlE  .= (view(in_rad.E_direct , iWLE, 1) .+ view(in_rad.E_diffuse, iWLE, 1)) .* dWL_iWLE;
    if photon
        sf_con.tmp_dwl_iWlE .*= WLE .* _FAC(FT);
    end;

    # calculate the SIF spectra for direct light
    # sf_con.M⁺ .= (Mb .+ Mf) ./ 2;
    # sf_con.M⁻ .= (Mb .- Mf) ./ 2;
    # mul!(sf_con.M⁺_sun, sf_con.M⁺, sf_con.tmp_dwl_iWlE);
    # mul!(sf_con.M⁻_sun, sf_con.M⁻, sf_con.tmp_dwl_iWlE);
    mul!(sf_con.M⁺_sun, Mb, sf_con.tmp_dwl_iWlE);
    mul!(sf_con.M⁻_sun, Mf, sf_con.tmp_dwl_iWlE);
    if photon
        sf_con.M⁺_sun ./= WLF .* _FAC(FT);
        sf_con.M⁻_sun ./= WLF .* _FAC(FT);
    end;

    # divide by pi to account for scattering
    sf_con.M⁻_sun .*= fqe / pi;
    sf_con.M⁺_sun .*= fqe / pi;

    return nothing
end
