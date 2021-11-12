###############################################################################
#
# leaf level simulations of PAR and APAR
#
###############################################################################
"""
    leaf_fluxes(leaf::LeafBios{FT},
                in_rad::IncomingRadiation{FT},
                wls::WaveLengths{FT},
                rt_con::RTCache{FT}
    ) where {FT<:AbstractFloat}

Return leaf PAR and APAR, given
- `leaf` [`LeafBios`](@ref) type struct
- `in_rad` [`IncomingRadiation`](@ref) type struct
- `wls` [`WaveLengths`](@ref) type struct
- `rt_con` [`RTCache`](@ref) type cache

Note that `in_rad` assumes direct light with zenith angle of 0, and a zenith
    angle correction needs to be made before passing it to this function.
"""
function leaf_fluxes(
            leaf::LeafBios{FT},
            in_rad::IncomingRadiation{FT},
            wls::WaveLengths{FT},
            rt_con::RTCache{FT}
) where {FT<:AbstractFloat}
    # unpack variables
    @unpack dWL_iPAR, iPAR, WL_iPAR = wls;
    cf_con = rt_con.cf_con;
    cf_con.kChlrel .= view(leaf.kChlrel, iPAR);

    # PAR energy from direct light
    cf_con.E_iPAR .= view(in_rad.E_direct, iPAR) .* view(leaf.α_SW, iPAR);
    e2phot!(WL_iPAR, cf_con.E_iPAR, cf_con.PAR_dir);

    # PAR energy from diffuse light
    cf_con.E_iPAR .= view(in_rad.E_diffuse, iPAR) .* view(leaf.α_SW, iPAR);
    e2phot!(WL_iPAR, cf_con.E_iPAR, cf_con.PAR_diff);

    # absorbed PAR energy from direct and diffuse light
    cf_con.PAR_diffCab .= cf_con.kChlrel .* cf_con.PAR_diff;
    cf_con.PAR_dirCab  .= cf_con.kChlrel .* cf_con.PAR_dir;

    # toral PAR and APAR
    _dif    = numerical∫(cf_con.PAR_diff   , dWL_iPAR);
    _dir    = numerical∫(cf_con.PAR_dir    , dWL_iPAR);
    _difCab = numerical∫(cf_con.PAR_diffCab, dWL_iPAR);
    _dirCab = numerical∫(cf_con.PAR_dirCab , dWL_iPAR);

    # return PAR and APAR in μmol m⁻² s⁻¹
    return 1000 * (_dir + _dif), 1000 * (_dirCab + _difCab)
end
