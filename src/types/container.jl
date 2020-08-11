###############################################################################
#
# Container of values to speed up the calculations
#
###############################################################################
"""
    mutable struct AngleContainer{FT}

A struct for canopy radiation information

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct AngleContainer{FT}
    "cos leaf azimuth angles"
    cos_ttlo ::Array{FT,1}
    "cos leaf azimuth angles - psi"
    cos_philo::Array{FT,1}
    "cos of normal of upperside of leaf"
    cos_ttli ::Array{FT,1}
    "sine of normal of upperside of leaf"
    sin_ttli ::Array{FT,1}
    "voscatt result"
    vol_scatt::Array{FT,1}

    # containers for more angles to use in canopy_geometry!
    "_Cs = cos_ttli .* cos(tts); # [nli]"
    _Cs::Array{FT,1}
    "_Ss = sin_ttli .* sin(tts); # [nli]"
    _Ss::Array{FT,1}
    "_Co = cos_ttli .* cos(tto); # [nli]"
    _Co::Array{FT,1}
    "_So = sin_ttli .* sin(tto); # [nli]"
    _So::Array{FT,1}
    "_1s = ones(FT,1,length(lazitab));"
    _1s::Array{FT,2}
    "cds = _Cs * _1s .+ _Ss * cos_ttlo';  # [nli, nlazi]"
    cds::Array{FT,2}
    "cdo = _Co * _1s .+ _So * cos_philo'; # [nli, nlazi]"
    cdo::Array{FT,2}
    "2D array to speed up; # [nli, nlazi]"
    _2d::Array{FT,2}
end
