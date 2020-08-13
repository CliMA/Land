###############################################################################
#
# Incoming radiation information
#
###############################################################################
"""
    struct IncomingRadiation{FT}

Incoming radiation information.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct IncomingRadiation{FT}
    "Wavelength `[nm]`"
    wl       ::Array{FT,1}
    "Direct incoming radiation `[mW m⁻² nm⁻¹]`"
    E_direct ::Array{FT,1}
    "Diffuse incoming radiation `[mW m⁻² nm⁻¹]`"
    E_diffuse::Array{FT,1}
end
