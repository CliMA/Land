"""
    canopy_fluxes!(can::Canopy4RT{FT}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, in_rad::IncomingRadiation{FT}, soil_opt::SoilOpticals{FT}, leaf_array::Array{LeafBios{FT},1}, wl_set::WaveLengths{FT}) where {FT<:AbstractFloat}

Computes a variety of integrated fluxes from the spectrally resolved computations in the short-wave Canopy RT (e.g. absorbed soil radiation, absorbed direct and diffuse PAR by layer (and angles for direct), net direct and diffuse energy balance per layer)
- `can` A [`Canopy4RT`](@ref) struct
- `can_opt` A [`CanopyOptiArray`](@ref) struct
- `can_rad` A [`CanopyRadiation`](@ref) struct
- `in_rad` An [`IncomingRadiationArray`](@ref) struct
- `soil_opt` A [`SoilOpti`](@ref) type struct for soil optical properties
- `leaf_array` An array of [`LeafBioArray`](@ref) type struct (i.e. leaf optical properties can change with canopy height)
- `wl_set` An [`WLParaSetArray`](@ref) type struct
"""
function canopy_fluxes!(
            can::Canopy4RT{FT},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            in_rad::IncomingRadiation{FT},
            soil_opt::SoilOpticals{FT},
            leaf_array::Array{LeafBios{FT},1},
            wl_set::WaveLengths{FT}
) where {FT<:AbstractFloat}
    @unpack dwl, iPAR, wl = wl_set

    # convert unit from mW to W
    fac             = FT(1e-3)
    nl              = can.nLayer
    iLAI            = (can.LAI*can.Ω)/nl;
    rsoil           = soil_opt.albedo_SW;
    soil_emissivity = 1 .-rsoil # emissivity of soil (goes into structure later)


    # Compute some fluxes, can be done separately if needed (this is absolute fluxes now, for the entire soil):
    can_rad.RnSoil_diffuse = fac * fast∫(dwl, can_rad.E_down[:,end].*soil_emissivity);
    can_rad.RnSoil_direct  = fac *fast∫(dwl, can_opt.Es_[:,end].*soil_emissivity);
    can_rad.RnSoil         = can_rad.RnSoil_direct + can_rad.RnSoil_diffuse

    # Normalization factor for leaf direct PAR (weighted sum has to be 1 to conserve net SW direct)
    normi = 1/mean(can_opt.absfs'*can.lidf)
    lPs = (can_opt.Ps[1:nl]+can_opt.Ps[2:nl+1])/2
    @inbounds for j = 1:nl
        if length(leaf_array)>1
            	kChlrel = leaf_array[j].kChlrel[iPAR]
        else
        kChlrel  = leaf_array[1].kChlrel[iPAR]
        end
        PAR_diff = (fac/iLAI)*e2phot(wl[iPAR], can_rad.netSW_shade[iPAR,j])
        # Direct PAR is normalized by layer Ps value:
        PAR_dir  = (fac/iLAI/lPs[j])*e2phot(wl[iPAR], can_rad.netSW_sunlit[iPAR,j])+PAR_diff

        PAR_diffCab = kChlrel.*PAR_diff
        #println(leaf.kChlrel[iPAR])
        # Direct PAR is normalized by layer Ps value:
        PAR_dirCab  = kChlrel.*PAR_dir;

        # Absorbed PAR per leaf for shaded:
        can_rad.absPAR_shade[j]      = fast∫(dwl[iPAR],PAR_diff);
        # for sunlit (all angles, normalized)
        #@show normi*can_opt.absfs*fast∫(dwl[iPAR],PAR_dir)
        can_rad.absPAR_sun[:,:,j]    = normi*can_opt.absfs*fast∫(dwl[iPAR],PAR_dir);
        # Same for true absorbed PAR by CabCar only  (needs to be taken into account in photosynthesis calculation, i.e. already accounts for leaf green fAPAR!)
        can_rad.absPAR_shadeCab[j]   = fast∫(dwl[iPAR],PAR_diffCab);
        can_rad.absPAR_sunCab[:,:,j] = normi*can_opt.absfs*fast∫(dwl[iPAR],PAR_dirCab);
    #println(PAR_dir[1], " ", can_rad.Pnh[j])
    end
    #@time fast∫(dwl[iPAR], in_rad.E_direct[iPAR])
    can_rad.incomingPAR_direct     = fac * fast∫(dwl[iPAR], e2phot(wl[iPAR],(in_rad.E_direct[iPAR])));
    can_rad.incomingPAR_diffuse    = fac * fast∫(dwl[iPAR], e2phot(wl[iPAR],(in_rad.E_diffuse[iPAR])));
    can_rad.incomingPAR            = can_rad.incomingPAR_diffuse*can_rad.incomingPAR_direct
    @inbounds for i=1:nl
        can_rad.intNetSW_shade[i]  = (fac/iLAI) *fast∫(dwl, can_rad.netSW_shade[:,i]);
        can_rad.intNetSW_sunlit[i] = (fac/iLAI/lPs[i]) * fast∫(dwl, can_rad.netSW_sunlit[:,i]) + can_rad.intNetSW_shade[i];
    end

    return nothing
end
