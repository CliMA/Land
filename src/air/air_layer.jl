#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-14: rename P_A to P_AIR, P_O2 to P_O₂
#     2022-Jan-14: add p_CO₂ to the structure
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
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jan-14: Move the structure from Photosynthesis.jl, only P_A and P_O2 for now
#     2022-Jan-24: add p_CO₂ to the constructors#
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
AirLayer{FT}() where {FT<:AbstractFloat} = AirLayer{FT}(P_ATM(FT), P_ATM(FT) * 0.209, 40);
