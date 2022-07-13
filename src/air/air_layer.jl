#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-14: rename P_A to P_AIR, P_O2 to P_O₂
#     2022-Jan-14: add p_CO₂ to the structure
#     2022-Jan-24: fix documentation
#     2022-Mar-09: add t, p_H₂O, p_H₂O_sat, rh, and wind fields
#     2022-Apr-19: move p_H₂O and wind to prognostic fields
#     2022-Jul-13: use @kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct AirLayer{FT<:AbstractFloat}
    # parameters that do not change with time
    "Atmospheric pressure `[Pa]`"
    P_AIR::FT = P_ATM()
    "O₂ partial pressure `[Pa]`"
    P_O₂::FT = P_ATM() * 0.209

    # prognostic variables that change with time
    "CO₂ partial pressure `[Pa]`"
    p_CO₂::FT = 40
    "H₂O partial pressure `[Pa]`"
    p_H₂O::FT = 1500
    "Temperature"
    t::FT = T_25()
    "Wind speed `[m s⁻¹]`"
    wind::FT = 1

    # diagnodtic variables that change with time
    "Saturated H₂O partial pressure `[Pa]`"
    p_H₂O_sat::FT = saturation_vapor_pressure(t)
    "relative humidity"
    rh::FT = p_H₂O / p_H₂O_sat
end
