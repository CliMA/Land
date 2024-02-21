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
    opti_file::String = LAND_2021

    # these should match the WL in wls
    "Direct incoming radiation `[mW m⁻² nm⁻¹]`"
    E_direct::Vector{FT} = read_nc(opti_file, "E_DIR")
    "Diffuse incoming radiation `[mW m⁻² nm⁻¹]`"
    E_diffuse::Vector{FT} = read_nc(opti_file, "E_DIFF")
end

IncomingRadiation(wls::WaveLengths{FT}) where {FT} = IncomingRadiation{FT}(opti_file = wls.opti_file)
