#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jun-02: migrate from CanopyLayers
#     2022-Jun-02: rename SolarAngles to SunSensorGeometry
#     2022-Jun-02: add extra fields HAA, HSA, saa, and vaa, and remove field raa (relative azimuth angle, will be computed based on saa and vaa)
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Structure that stores sun sensor geometry information.

# Fields

$(TYPEDFIELDS)

"""
mutable struct SunSensorGeometry{FT<:AbstractFloat}
    # parameters that do not change with time
    "Hill facing azimuth angle `[°]`, 0 for North, 180 for south"
    HAA::FT
    "Hill slope angle `[°]`"
    HSA::FT

    # parameters that do not change with time
    "Solar azimuth angle `[°]`, a function of time"
    saa::FT
    "Solar zenith angle `[°]`, a function of lat and time"
    sza::FT
    "Viewing azimuth angle `[°]`"
    vaa::FT
    "Viewing zenith angle `[°]`"
    vza::FT
end


#######################################################################################################################################################################################################
#
# Changes to this constructor
# General
#     2022-Jun-02: add constructor
#
#######################################################################################################################################################################################################
"""

    SunSensorGeometry{FT}(; haa::Number = 0, hsa::Number = 0, saa::Number = 180, sza::Number = 30, vaa::Number = 180, vza::Number = 0) where {FT<:AbstractFloat}

Constructor for [`SunSensorGeometry`](@ref), given
- `haa` Hill facing azimuth angle
- `hsa` Hill slope angle
- `saa` Solar azimuth angle
- `sza` Solar zenith angle
- `vaa` Viewing azimuth angle
- `vza` Viewing zenith angle
"""
SunSensorGeometry{FT}(; haa::Number = 0, hsa::Number = 0, saa::Number = 180, sza::Number = 30, vaa::Number = 180, vza::Number = 0) where {FT<:AbstractFloat} = (
    return SunSensorGeometry{FT}(
                haa,    # HAA
                hsa,    # HSA
                saa,    # saa
                sza,    # sza
                vaa,    # vaa
                vza     # vza
    )
);
