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
    cos_ttlo ::Array{FT,1} = cosd.(lazitab)
    "cos leaf azimuth angles - psi"
    cos_philo::Array{FT,1} = cosd.(lazitab .- psi)
    "cos of normal of upperside of leaf"
    cos_ttli ::Array{FT,1} = cosd.(litab)
    "sine of normal of upperside of leaf"
    sin_ttli ::Array{FT,1} = sind.(litab)
end
