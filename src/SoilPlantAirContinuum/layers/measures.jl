"""
    GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the GPP of the SPAC per ground area
"""
function GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _gpp::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _gpp += numerical∫(_iPS.Ag, _iPS.LAIx) * _iPS.LA;
    end;

    return _gpp / spac.ga
end


"""
    CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the canopy NPP of the SPAC per ground area
"""
function CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _cnpp::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _cnpp += numerical∫(_iPS.An, _iPS.LAIx) * _iPS.LA;
    end;

    return _cnpp / spac.ga
end


"""
    T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the transpiration of the SPAC per ground area
"""
function T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _t_veg::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iEN = spac.envirs[_i_can];
        _iPS = spac.plant_ps[_i_can];
        _t_veg += numerical∫(_iPS.g_lw, _iPS.LAIx) * (_iPS.p_sat - _iEN.p_H₂O) / _iEN.p_atm * _iPS.LA;
    end;

    return _t_veg / spac.ga
end


"""
    PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the cumulative PPAR of the SPAC per ground area
"""
function PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    _ppar::FT = 0;

    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can];
        _ppar += numerical∫(_iPS.APAR, _iPS.LAIx) * FT(1e-6) * spac.canopy_rt.LAI / spac.canopy_rt.nLayer;
    end;

    return _ppar
end
