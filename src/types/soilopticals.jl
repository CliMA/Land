###############################################################################
#
# Soil optical parameters
#
###############################################################################
"""
    mutable struct SoilOpticals{FT}

A struct of soil optical parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct SoilOpticals{FT}
    "Wavelength `[nm]`"
    WL        ::Array{FT,1}
    "Shortwave albedo"
    albedo_SW ::Array{FT,1}
    "Shortwave Emissivity"
    emsvty_SW ::Array{FT,1}
    "Longwave albedo"
    albedo_LW ::Array{FT,1}
    "Soil surface temperature `[K]`"
    soil_skinT::FT
end
