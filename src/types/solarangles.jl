###############################################################################
#
# Solar angle type
#
###############################################################################
"""
    struct SolarAngles{FT}

Struct for observation and solar angles

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SolarAngles{FT}
    "Hill zenith angle `[°]`, hill slope angle"
    hza::FT = 0
    "Hill azimuth angle `[°]`, 0 for North, 180 for south"
    haa::FT = 180
    "Solar azimuth angle `[°]`, a function of time"
    saa::FT = 0
    "Solar zenith angle `[°]`, a function of lat and time"
    sza::FT = 30
    "Viewing azimuth angle `[°]`"
    vaa::FT = 0
    "Viewing zenith angle `[°]`"
    vza::FT = 0
    "Relative azimuth angle `[°]`, difference between saa and vaa"
    raa::FT = 0
end
