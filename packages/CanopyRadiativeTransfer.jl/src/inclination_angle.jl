#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: generalize the function from CanopyLayers.dcum to lidf_cdf
#
#######################################################################################################################################################################################################
"""

    lidf_cdf(lidf::VerhoefLIDF{FT}, θ::FT) where {FT<:AbstractFloat}

Return the cumulative distribution frequency, given
- `lidf` `VerhoefLIDF` type algorithm
- `θ` Leaf inclination angle in `[°]`

"""
function lidf_cdf end

lidf_cdf(lidf::VerhoefLIDF{FT}, θ::FT) where {FT<:AbstractFloat} = (
    (; A, B) = lidf;

    if A >= 1
        return 1 - cosd(θ)
    end;

    # iterate to solve for CDF solution
    _θ = deg2rad(θ);
    _y::FT = 0;
    _x = 2 * _θ;
    _n = 0;
    _δx::FT = 1;
    while (abs(_δx) >= max(eps(FT), 1e-8)) && (_n < 50)
        _y = A * sin(_x) + B / 2 * sin(2*_x);
        _δx = (_y - _x) / 2 + _θ;
        _x += _δx;
        _n += 1;
    end;

    return 2 * (_y + _θ) / FT(π)
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jun-02: generalize the function from CanopyLayers.dladgen to inclination_angles!
#     2022-Jun-02: add method for VerhoefLIDF algorithm
#     2022-Jun-15: add support to BroadbandSLCanopy
#
#######################################################################################################################################################################################################
"""

    inclination_angles!(can::Union{BroadbandSLCanopy{FT}, HyperspectralMLCanopy{FT}}, lidf::VerhoefLIDF{FT}) where {FT<:AbstractFloat}

Update the frequency of leaf inclination angles, given
- `can` `HyperspectralMLCanopy` type multiple layer canopy
- `lidf` `VerhoefLIDF` type algorithm

"""
function inclination_angles! end

inclination_angles!(can::Union{BroadbandSLCanopy{FT}, HyperspectralMLCanopy{FT}}, lidf::VerhoefLIDF{FT}) where {FT<:AbstractFloat} = (
    (; Θ_INCL_BNDS) = can;

    for _i in eachindex(can.P_INCL)
        can.P_INCL[_i] = lidf_cdf(lidf, Θ_INCL_BNDS[_i,2]) - lidf_cdf(lidf, Θ_INCL_BNDS[_i,1]);
    end;

    return nothing
);
