###############################################################################
#
# Incoming radiation information
#
###############################################################################
"""
    mutable struct IncomingRadiation{FT}

Incoming radiation information

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct IncomingRadiation{FT}
    # these should match the WL in wls
    "Direct incoming radiation `[mW m⁻² nm⁻¹]`"
    E_direct ::Array{FT,1}
    "Diffuse incoming radiation `[mW m⁻² nm⁻¹]`"
    E_diffuse::Array{FT,1}
end
