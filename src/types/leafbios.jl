###############################################################################
#
# Leaf biological parameters
#
###############################################################################
"""
    struct LeafBios{FT}

A struct of leaf biological parameters

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct LeafBios{FT}
    # Array size information
    "Number of wave length"
    nWL ::Int = 10
    "Number of wave length for excitation"
    nWLe::Int = 10
    "Number of wave length for SIF"
    nWLf::Int = 10

    "Leaf structure parameter"
    N   ::FT = FT(1.4  )
    "Chlorophyll a+b content `[µg cm⁻²]`"
    Cab ::FT = FT(40.0 )
    "Carotenoid content `[µg cm⁻²]`"
    Car ::FT = FT(10.0 )
    "Anthocynanin content `[µg cm⁻²]`"
    Ant ::FT = FT(0.0  )
    "Senescent material fraction"
    Cs  ::FT = FT(0.0  )
    "Equivalent water thickness `[cm]`"
    Cw  ::FT = FT(0.009)
    "Dry matter content (dry leaf mass per unit area) `[g cm⁻²]`"
    Cm  ::FT = FT(0.012)
    "Fractionation between Zeaxanthin and Violaxanthin in Car (1=all Zeaxanthin) (-)"
    Cx  ::FT = FT(0.0  )
    "Broadband thermal reflectance (-)"
    ρ_LW::FT = FT(0.01 )
    "Broadband thermal transmission (-)"
    τ_LW::FT = FT(0.01 )
    "Leaf fluorescence efficiency (Fo standard)"
    fqe ::FT = FT(1.0  )

    "Shortwave leaf reflectance"
    ρ_SW       ::Array{FT,1} = zeros(FT, nWL)
    "Shortwave leaf transmission"
    τ_SW       ::Array{FT,1} = zeros(FT, nWL)
    "Relative absorbtion by Chlorophyll+Car"
    kChlrel    ::Array{FT,1} = zeros(FT, nWL)
    "Relative absorbtion by Chlorophyll"
    kChlrel_old::Array{FT,1} = zeros(FT, nWL)
    "Fluorescence excitation matrix backwards"
    Mb         ::Array{FT,2} = zeros(FT,(nWLf,nWLe))
    "Fluorescence excitation matrix forwards"
    Mf         ::Array{FT,2} = zeros(FT,(nWLf,nWLe))
    "Doubling adding layers"
    ndub      ::Int = 10
end




"""
    create_leaf_bios(FT, nWl::Int, nWle::Int, nWlf::Int)

Create a leaf biological parameters struct, given
- `FType` Floating number type
- `nWl` Number of wave length
- `nWle` Number of excitation wave length
- `nWlf` Number of fluorescence wave length

Returns a [`LeafBioMArray`](@ref) or [`LeafBioArray`](@ref) type struct.
"""
function create_leaf_bios(FT, nWL::Int, nWLe::Int, nWLf::Int)
    return LeafBios{FT}(nWL=nWL, nWLe=nWLe, nWLf=nWLf)
end
