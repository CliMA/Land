###############################################################################
#
# Update canopy fluxes
#
###############################################################################
"""
    canopy_fluxes!(can::Canopy4RT{FT}, can_opt::CanopyOpticals{FT}, can_rad::CanopyRads{FT}, in_rad::IncomingRadiation{FT}, soil_opt::SoilOpticals{FT}, leaf_array::Array{LeafBios{FT},1}, wl_set::WaveLengths{FT}, rt_con::RTContainer) where {FT<:AbstractFloat}

Computes a variety of integrated fluxes from the spectrally resolved
    computations in the short-wave Canopy RT (e.g. absorbed soil radiation,
    absorbed direct and diffuse PAR by layer (and angles for direct), net
    direct and diffuse energy balance per layer), given
- `can` A [`Canopy4RT`](@ref) struct
- `can_opt` A [`CanopyOpticals`](@ref) struct
- `can_rad` A [`CanopyRads`](@ref) struct
- `in_rad` An [`IncomingRadiation`](@ref) struct
- `soil_opt` A [`SoilOpticals`](@ref) type struct for soil optical properties
- `leaf_array` An array of [`LeafBios`](@ref) type struct (i.e. leaf optical
    properties can change with canopy height)
- `wl_set` An [`WaveLengths`](@ref) type struct
- `rt_con` [`RTContainer`](@ref) type container
"""
function canopy_fluxes!(
            can::Canopy4RT{FT},
            can_opt::CanopyOpticals{FT},
            can_rad::CanopyRads{FT},
            in_rad::IncomingRadiation{FT},
            soil_opt::SoilOpticals{FT},
            leaf_array::Array{LeafBios{FT},1},
            wl_set::WaveLengths{FT},
            rt_con::RTContainer{FT}
) where {FT<:AbstractFloat}
    # 1. unpack variables from structures
    @unpack LAI, nLayer, Ω = can;
    @unpack albedo_SW, emsvty_SW = soil_opt;
    @unpack dWL, iPAR, WL = wl_set;

    # 2. compute some useful variables
    fac  = FT(1e-3);
    iLAI = LAI * Ω / nLayer;

    # 3. Compute some fluxes, can be done separately if needed
    #    this is absolute fluxes now, for the entire soil
    last_ind_cr            = lastindex(can_rad.E_down,2);
    rt_con.abs_wave       .= view(can_rad.E_down, :, last_ind_cr) .* emsvty_SW;
    can_rad.RnSoil_diffuse = fac * fast∫!(rt_con.abs_wave, dWL);
    rt_con.abs_wave       .= view(can_opt.Es_, :, last_ind_cr) .* emsvty_SW;
    can_rad.RnSoil_direct  = fac * fast∫!(rt_con.abs_wave, dWL);
    can_rad.RnSoil         = can_rad.RnSoil_direct + can_rad.RnSoil_diffuse;

    # 4. Normalization factor for leaf direct PAR
    #    weighted sum has to be 1 to conserve net SW direct
    #    Direct PAR is normalized by layer Ps value
    mul!(rt_con.absfs_lidf, can_opt.absfs', can.lidf);
    normi       = 1 / mean(rt_con.absfs_lidf);
    rt_con.lPs .= (view(can_opt.Ps, 1:nLayer  ) .+
                  view(can_opt.Ps, 2:nLayer+1)) ./ 2;
    @unpack dλ_iPAR, lPs, λ_iPAR = rt_con;

    @inbounds for j in 1:nLayer
        if length(leaf_array)>1
            rt_con.kChlrel .= view(leaf_array[j].kChlrel, iPAR);
        else
            rt_con.kChlrel .= view(leaf_array[1].kChlrel, iPAR);
        end
        # for diffuse PAR
        rt_con.E_iPAR .= view(can_rad.netSW_shade, iPAR, j);
        e2phot!(λ_iPAR, rt_con.E_iPAR, rt_con.PAR_diff);
        rt_con.PAR_diff .*= fac / iLAI;
        # for direct PAR
        rt_con.E_iPAR .= view(can_rad.netSW_sunlit, iPAR, j);
        e2phot!(λ_iPAR, rt_con.E_iPAR, rt_con.PAR_dir);
        rt_con.PAR_dir .*= fac / iLAI;
        rt_con.PAR_dir .+= rt_con.PAR_diff;
        # for leaf absorbed
        rt_con.PAR_diffCab .= rt_con.kChlrel .* rt_con.PAR_diff;
        rt_con.PAR_dirCab  .= rt_con.kChlrel .* rt_con.PAR_dir;

        # Absorbed PAR per leaf for shaded, PAR_diff and PAR_dir changed!
        _dif = fast∫!(rt_con.PAR_diff, dλ_iPAR);
        _dir = fast∫!(rt_con.PAR_dir, dλ_iPAR) * normi;
        can_rad.absPAR_shade[j] = _dif;
        can_rad.absPAR_sun[:,:,j] .= can_opt.absfs .* _dir;

        # absorbed PAR for photosynthesis
        _difCab = fast∫!(rt_con.PAR_diffCab, dλ_iPAR);
        _dirCab = fast∫!(rt_con.PAR_dirCab , dλ_iPAR) * normi;
        can_rad.absPAR_shadeCab[j] = _difCab;
        can_rad.absPAR_sunCab[:,:,j] .= can_opt.absfs .* _dirCab;
    end

    # 5. Total PAR
    # TODO considering remove this part, if we are not using it
    #    re-use the PAR_dir and PAR_diff in the rt_con
    rt_con.E_iPAR .= view(in_rad.E_direct, iPAR);
    e2phot!(λ_iPAR, rt_con.E_iPAR, rt_con.PAR_dir);
    rt_con.E_iPAR .= view(in_rad.E_diffuse, iPAR);
    e2phot!(λ_iPAR, rt_con.E_iPAR, rt_con.PAR_diff);
    can_rad.incomingPAR_direct  = fac * fast∫!(rt_con.PAR_dir , dλ_iPAR);
    can_rad.incomingPAR_diffuse = fac * fast∫!(rt_con.PAR_diff, dλ_iPAR);
    can_rad.incomingPAR         = can_rad.incomingPAR_diffuse +
                                  can_rad.incomingPAR_direct;
    @inbounds for i in 1:nLayer
        rt_con.E_all .= view(can_rad.netSW_shade, :, i);
        can_rad.intNetSW_shade[i]  = fast∫!(rt_con.E_all, dWL) * fac / iLAI;
        rt_con.E_all .= view(can_rad.netSW_sunlit, :, i);
        can_rad.intNetSW_sunlit[i] = fast∫!(rt_con.E_all, dWL) *
                                     fac / iLAI / lPs[i] +
                                     can_rad.intNetSW_shade[i];
    end

    return nothing
end
