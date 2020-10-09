###############################################################################
#
# Solar angle type
#
###############################################################################
"""
    struct SolarAngles{FT}

Struct for observation and solar angles

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct SolarAngles{FT}
    "Solar Zenith Angle `[degree]`"
    tts::FT = FT(30.0)
    "Viewing Zenith Angle in `[degree]`"
    tto::FT = FT(0.0 )
    "relative azimuth in `[degree]`"
    psi::FT = FT(0.0 )
end
