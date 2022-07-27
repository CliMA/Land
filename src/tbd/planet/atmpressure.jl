###############################################################################
#
# Atmospheric pressure and partial pressure
#
###############################################################################
GM2RT25(FT) = GRAVITY(FT) * FT(0.02896968) / RT_25(FT);

"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure_ratio(
            h::FT
) where {FT<:AbstractFloat}
    return exp( -h * GM2RT25(FT) );
end




"""
    atmospheric_pressure(h::FT)

Calculate the atmospheric pressure, given
- `h` elevation in `[m]`
"""
function atmospheric_pressure(
            h::FT
) where {FT<:AbstractFloat}
    return P_ATM(FT) * atmospheric_pressure_ratio(h);
end
