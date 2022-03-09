#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-14: rename P_A to P_AIR, P_O2 to P_O₂
#     2022-Jan-14: add p_CO₂ to the structure
#     2022-Jan-24: fix documentation
#     2022-Mar-09: add t, p_H₂O, p_H₂O_sat, rh, and wind fields
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores air layer information

# Fields

$(TYPEDFIELDS)

"""
mutable struct AirLayer{FT<:AbstractFloat}
    # parameters that do not change with time
    "Atmospheric pressure `[Pa]`"
    P_AIR::FT
    "O₂ partial pressure `[Pa]`"
    P_O₂::FT

    # prognostic variables that change with time
    "CO₂ partial pressure `[Pa]`"
    p_CO₂::FT
    "Temperature"
    t::FT

    # diagnodtic variables that change with time
    "H₂O partial pressure `[Pa]`"
    p_H₂O::FT
    "Saturated H₂O partial pressure `[Pa]`"
    p_H₂O_sat::FT
    "relative humidity"
    rh::FT
    "Wind speed `[m s⁻¹]`"
    wind::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-24: add p_CO₂ to the constructors
#     2022-Mar-09: add t, p_H₂O, p_H₂O_sat, rh, and wind fields
#
#######################################################################################################################################################################################################
"""

    AirLayer{FT}() where {FT<:AbstractFloat}

Constructor for AirLayer

---
# Examples
```julia
air = AirLayer{Float64}();
```
"""
AirLayer{FT}() where {FT<:AbstractFloat} = (
    return AirLayer{FT}(
                P_ATM(FT),                                  # P_AIR
                P_ATM(FT) * 0.209,                          # P_O₂
                40,                                         # p_CO₂
                T_25(FT),                                   # t
                1500,                                       # p_H₂O
                saturation_vapor_pressure(T_25(FT)),        # p_H₂O_sat
                1500 / saturation_vapor_pressure(T_25(FT)), # rh
                1)                                          # wind
);
