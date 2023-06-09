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
    nr::Vector{FT} = zeros(FT, 2)
    Km::Vector{FT} = zeros(FT, 2)
    Kab::Vector{FT} = zeros(FT, 2)
    Kant::Vector{FT} = zeros(FT, 2)
    Kcar::Vector{FT} = zeros(FT, 2)
    Kw::Vector{FT} = zeros(FT, 2)
    KBrown::Vector{FT} = zeros(FT, 2)
    phi::Vector{FT} = zeros(FT, 2)
    KcaV::Vector{FT} = zeros(FT, 2)
    KcaZ::Vector{FT} = zeros(FT, 2)
    "Wave length `[nm]`, same as `WL` in [`WaveLengths`](@ref)`"
    lambda::Vector{FT} = zeros(FT, 2)
end
