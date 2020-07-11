###############################################################################
#
# Atmospheric pressure and partial pressure
#
###############################################################################
const GM2RT25 = GRAVITY * 0.02896968 / RK_25;

"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure_ratio(
            h::FT
) where {FT<:AbstractFloat}
    return exp( -h * FT(GM2RT25) );
end




"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure(
            h::FT
) where {FT<:AbstractFloat}
    return FT(P_ATM) * atmospheric_pressure_ratio(h);
end




"""
    ppm_to_Pa(h::FT)

Convert ppm to Pa, given
- `h` elevation in `[m]`
"""
function ppm_to_Pa(
            h::FT
) where {FT<:AbstractFloat}
    return atmospheric_pressure(h) * FT(1e-6)
end

