#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jul-14: add Meteorology struct to store meteorological data
#     2022-Jul-14: add field t_precip
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores meteorological information

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Meteorology{FT<:AbstractFloat}
    # prognostic variables that change with time
    "Precipitation in form of rain (before interception) `[mol m⁻²]`"
    rain::FT = 0
    "Precipitation in form of snow (before interception) `[mol m⁻²]`"
    snow::FT = 0
    "Precipitation temperature `[K]`"
    t_precip::FT = T_25()
end
