###############################################################################
#
# Reflectance related
#
###############################################################################
"""
    BLUE(can_rad::CanopyRads{FT},
         wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the BLUE @ 469 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function BLUE(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return REF_WL(can_rad, wls, FT(469))
end




"""
    EVI(can_rad::CanopyRads{FT},
        wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the EVI, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function EVI(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    _BLUE = BLUE(can_rad, wls);
    _NIR  = NIR(can_rad, wls);
    _RED  = RED(can_rad, wls);

    return FT(2.5) * (_NIR - _RED) / (_NIR + 6 * _RED - FT(7.5) * _BLUE + 1)
end




"""
    EVI2(can_rad::CanopyRads{FT},
        wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the EVI2, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function EVI2(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    _BLUE = BLUE(can_rad, wls);
    _NIR  = NIR(can_rad, wls);
    _RED  = RED(can_rad, wls);

    return FT(2.5) * (_NIR - _RED) / (_NIR + FT(2.4) * _RED + 1)
end




"""
    LSWI(can_rad::CanopyRads{FT},
         wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the LSWI, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function LSWI(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    _BLUE = BLUE(can_rad, wls);
    _NIR  = NIR(can_rad, wls);
    _RED  = RED(can_rad, wls);
    _SWIR = SWIR(can_rad, wls);


    return (_NIR - _SWIR) / (_NIR + _SWIR)
end




"""
    NDVI(can_rad::CanopyRads{FT},
         wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the NDVI, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function NDVI(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    _NIR = NIR(can_rad, wls);
    _RED = RED(can_rad, wls);

    return (_NIR - _RED) / (_NIR + _RED)
end




"""
    NIR(can_rad::CanopyRads{FT},
        wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the NIR @ 858.5 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function NIR(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return REF_WL(can_rad, wls, FT(858.5))
end




"""
    NIRv(can_rad::CanopyRads{FT},
         wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the NIRv, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function NIRv(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    _NIR = NIR(can_rad, wls);
    _RED = RED(can_rad, wls);

    return (_NIR - _RED) / (_NIR + _RED) * _NIR
end




"""
    NIRvES(can_rad::CanopyRads{FT},
           wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the NIRv from the energy spectrum, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function NIRvES(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    _NIR = SPECTRUM_WL(can_rad.Lo, wls, FT(858.5));
    _RED = SPECTRUM_WL(can_rad.Lo, wls, FT(645));

    return (_NIR - _RED) / (_NIR + _RED) * _NIR
end




"""
    RED(can_rad::CanopyRads{FT},
        wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the RED @ 645 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function RED(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return REF_WL(can_rad, wls, FT(645))
end




"""
    REF_WL(wls::WaveLengths{FT},
           can_rad::CanopyRads{FT}
           wls::WaveLengths{FT},
           twl::FT
    ) where {FT<:AbstractFloat}

Return the Reflectance, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `twl` Target wave length in nm
"""
function REF_WL(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT},
            twl::FT
) where {FT<:AbstractFloat}
    @unpack WL,nWL = wls;
    @unpack alb_obs = can_rad;

    # find the index where twl nm is
    ind = 0;
    for i in 1:(nWL-1)
        if WL[i] <= twl <= WL[i+1]
            ind = i;
            break;
        end
    end

    # warning if ind == 0
    if ind==0
        @warn "target wave length out of bounds, please check the set up!";
        return FT(NaN)
    else
        return (WL[ind+1] - twl) / (WL[ind+1] - WL[ind]) * alb_obs[ind] +
               (twl - WL[ind]) / (WL[ind+1] - WL[ind]) * alb_obs[ind+1]
    end
end




"""
    SPECTRUM_WL(wls::WaveLengths{FT},
           can_rad::CanopyRads{FT}
           wls::WaveLengths{FT},
           twl::FT
    ) where {FT<:AbstractFloat}

Return the Reflectance, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `twl` Target wave length in nm
"""
function SPECTRUM_WL(
            spectrum::Vector{FT},
            wls::WaveLengths{FT},
            twl::FT
) where {FT<:AbstractFloat}
    @unpack WL,nWL = wls;

    # find the index where twl nm is
    ind = 0;
    for i in 1:(nWL-1)
        if WL[i] <= twl <= WL[i+1]
            ind = i;
            break;
        end
    end

    # warning if ind == 0
    if ind==0
        @warn "target wave length out of bounds, please check the set up!";
        return FT(NaN)
    else
        return (WL[ind+1] - twl) / (WL[ind+1] - WL[ind]) * spectrum[ind] +
               (twl - WL[ind]) / (WL[ind+1] - WL[ind]) * spectrum[ind+1]
    end
end




"""
    SWIR(can_rad::CanopyRads{FT},
        wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the SWIR @ 2130 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function SWIR(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return REF_WL(can_rad, wls, FT(2130))
end








###############################################################################
#
# SIF related
#
###############################################################################
"""
    SIF_WL(wls::WaveLengths{FT},
           can_rad::CanopyRads{FT}
           wls::WaveLengths{FT},
           twl::FT
    ) where {FT<:AbstractFloat}

Return the SIF, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `twl` Target SIF wave length in nm
"""
function SIF_WL(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT},
            twl::FT
) where {FT<:AbstractFloat}
    @unpack WLF,nWLF = wls;
    @unpack SIF_obs = can_rad;

    # find the index where twl nm is
    ind = 0;
    for i in 1:(nWLF-1)
        if WLF[i] <= twl <= WLF[i+1]
            ind = i;
            break;
        end
    end

    # warning if ind == 0
    if ind==0
        @warn "target wave length out of bounds, please check the set up!";
        return FT(NaN)
    else
        return (WLF[ind+1] - twl) / (WLF[ind+1] - WLF[ind]) * SIF_obs[ind] +
               (twl - WLF[ind]) / (WLF[ind+1] - WLF[ind]) * SIF_obs[ind+1]
    end
end




"""
    SIF_740(can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the SIF @ 740 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function SIF_740(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return SIF_WL(can_rad, wls, FT(740))
end




"""
    SIF_757(can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the SIF @ 757 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function SIF_757(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return SIF_WL(can_rad, wls, FT(757))
end




"""
    SIF_771(can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
    ) where {FT<:AbstractFloat}

Return the SIF @ 771 nm, given
- `can_rad` [`CanopyRads`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
"""
function SIF_771(
            can_rad::CanopyRads{FT},
            wls::WaveLengths{FT}
) where {FT<:AbstractFloat}
    return SIF_WL(can_rad, wls, FT(771))
end
