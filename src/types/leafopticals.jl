###############################################################################
#
# Leaf optical parameters
#
###############################################################################
"""
    mutable struct LeafOpticals{FT}

Struct for leaf optical properties

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct LeafOpticals{FT}
    # TODO Add explanations to each field
    nr    ::Array{FT,1} = zeros(FT, 2)
    Km    ::Array{FT,1} = zeros(FT, 2)
    Kab   ::Array{FT,1} = zeros(FT, 2)
    Kant  ::Array{FT,1} = zeros(FT, 2)
    Kcar  ::Array{FT,1} = zeros(FT, 2)
    Kw    ::Array{FT,1} = zeros(FT, 2)
    KBrown::Array{FT,1} = zeros(FT, 2)
    phi   ::Array{FT,1} = zeros(FT, 2)
    KcaV  ::Array{FT,1} = zeros(FT, 2)
    KcaZ  ::Array{FT,1} = zeros(FT, 2)
    "Wave length `[nm]`, same as `WL` in [`WaveLengths`](@ref)`"
    lambda::Array{FT,1} = zeros(FT, 2)
end
