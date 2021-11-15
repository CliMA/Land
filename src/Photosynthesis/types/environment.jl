###############################################################################
#
# Environmental conditions container
#
###############################################################################
"""
    mutable struct AirLayer{FT}

Struct to store environmental conditions in each air layer corresponds to one
    canopy layer.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    "Air temperature `[K]`"
    t_air::FT = 298.15

    # pressures
    "Atmospheric CO₂ partial pressure `[Pa]`"
    p_a::FT = FT(41.0)
    "Atmospheric pressure `[Pa]`"
    p_atm::FT = FT(101325.0)
    "Atmospheric vapor pressure `[Pa]`"
    p_H₂O::FT = FT(1500.0)
    "Atmospheric O₂ partial pressure `[Pa]`"
    p_O₂::FT = FT(101325.0 * 0.209)
    "Saturation vapor pressure `[Pa]`"
    p_sat::FT = saturation_vapor_pressure(t_air)
    "Relatiev humidity"
    RH::FT = p_H₂O / p_sat
    "Vapor pressure deficit `[Pa]`"
    vpd::FT = p_sat - p_H₂O

    # wind speed
    "Wind speed `[m s⁻¹]`"
    wind::FT = FT(2)
end
