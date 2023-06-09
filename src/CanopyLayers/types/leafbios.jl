###############################################################################
#
# Leaf biological parameters
#
###############################################################################
"""
    mutable struct LeafBios{FT}

A struct of leaf biological parameters

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct LeafBios{FT}
    # local storage of dimension information
    "Number of wave length"
    nWL::Int = 10
    "Number of wave length for excitation"
    nWLE::Int = 10
    "Number of wave length for SIF"
    nWLF::Int = 10

    "Leaf structure parameter"
    N::FT = 1.4
    "Chlorophyll a+b content `[µg cm⁻²]`"
    Cab::FT = 40
    "Carotenoid content `[µg cm⁻²]`"
    Car::FT = 10
    "Anthocynanin content `[µg cm⁻²]`"
    Ant::FT = 0
    "Senescent material fraction"
    Cs::FT = 0
    "Equivalent water thickness `[cm]`"
    Cw::FT = 0.009
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    Cm::FT = 0.012
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin) (-)"
    Cx::FT = 0
    "Leaf fluorescence efficiency (Fo standard)"
    fqe::FT = 1
    "Broadband thermal reflectance (-)"
    ρ_LW::FT = 0.01
    "Broadband thermal transmission (-)"
    τ_LW::FT = 0.01

    "Shortwave leaf reflectance"
    ρ_SW::Vector{FT} = zeros(FT, nWL)
    "Shortwave leaf transmission"
    τ_SW::Vector{FT} = zeros(FT, nWL)
    "Shortwave absorption"
    α_SW::Vector{FT}  = zeros(FT, nWL)
    "Relative absorbtion by Chlorophyll+Car"
    kChlrel::Vector{FT} = zeros(FT, nWL)
    "Relative absorbtion by Chlorophyll"
    kChlrel_old::Vector{FT} = zeros(FT, nWL)
    "Fluorescence excitation matrix backwards"
    Mb::Matrix{FT} = zeros(FT,(nWLF,nWLE))
    "Fluorescence excitation matrix forwards"
    Mf::Matrix{FT} = zeros(FT,(nWLF,nWLE))
    "Doubling adding layers"
    ndub::Int = 10
end

LeafBios{FT}(rt_dim::RTDimensions) where {FT} = (
    (; nWL, nWLE, nWLF) = rt_dim;

    return LeafBios{FT}(nWL = nWL, nWLE = nWLE, nWLF = nWLF)
);
