#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename SolarAngles to SunSensorGeometry
#     2022-Jun-02: add extra fields HAA, HSA, saa, and vaa, and remove field raa (relative azimuth angle, will be computed based on saa and vaa)
#     2022-Jul-20: use kwdef for the constructor
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores sun sensor geometry information.

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct SunSensorGeometry{FT<:AbstractFloat}
    # General site information
    "Hill facing azimuth angle `[°]`, 0 for North, 180 for south"
    HAA::FT = 0
    "Hill slope angle `[°]`"
    HSA::FT = 0

    # Prognostic variables
    "Solar azimuth angle `[°]`, a function of time"
    saa::FT = 180
    "Solar zenith angle `[°]`, a function of lat and time"
    sza::FT = 30
    "Viewing azimuth angle `[°]`"
    vaa::FT = 180
    "Viewing zenith angle `[°]`"
    vza::FT = 0
end
