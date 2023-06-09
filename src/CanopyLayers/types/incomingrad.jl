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
    E_direct::Vector{FT}
    "Diffuse incoming radiation `[mW m⁻² nm⁻¹]`"
    E_diffuse::Vector{FT}
end
