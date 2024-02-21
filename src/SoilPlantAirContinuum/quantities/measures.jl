"""
    T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the transpiration of the SPAC per ground area
"""
function T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    t_veg::FT = 0;

    for i_can in 1:spac.n_canopy
        iEN = spac.envirs[i_can];
        iPS = spac.plant_ps[i_can];
        t_veg += numerical∫(iPS.g_lw, iPS.LAIx) * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm * iPS.LA;
    end;

    return t_veg / spac.ga
end


"""
    PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the cumulative PPAR of the SPAC per ground area
"""
function PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    ppar::FT = 0;

    for iPS in spac.plant_ps
        ppar += numerical∫(iPS.APAR, iPS.LAIx) * FT(1e-6) * spac.canopy_rt.LAI / spac.canopy_rt.nLayer;
    end;

    return ppar
end
