###############################################################################
#
# Canopy4RT structure
#
###############################################################################
"""
    mutable struct Canopy4RT

A canopy struct for the radiation transfer module

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Canopy4RT{FT<:AbstractFloat}
    # local storage of dimension information
    "Number of canopy layers"
    nLayer::Int = 20

    # canopy information
    "Leaf Area Index"
    LAI::FT = 3
    "Clumping factor"
    Ω::FT = 1
    "Structure factor a"
    clump_a::FT = 1
    "Structure factor b"
    clump_b::FT = 0
    "Leaf width"
    leaf_width::FT = 0.1
    "Vegetation height"
    hc::FT = 2
    "Leaf Inclination"
    LIDFa::FT = FT(0)
    "Variation in leaf inclination"
    LIDFb::FT = FT(0)
    "HotSpot parameter (still need to check!)"
    hot::FT = 0.05

    # tree/canopy/leaf traits
    "Canopy height `[m]`"
    height::FT = 20
    "Canopy roughness `[m]`"
    z0m::FT = 1
    "Tree roughtnes `[m]`"
    z0h::FT = -999
    "Canopy displacement height `[m]`"
    d::FT = -999
    "m/sqrt(s) turbulent transfer coefficient"
    Cd::FT = 0.01

    # Some more derived parameters:
    "List of mean inclination angles `[°]`"
    litab::Vector{FT} = collect(FT,5:10:85)
    "List of inclination angle boundaries `[°]`"
    litab_bnd::Matrix{FT} = [collect(0:10:80) collect(FT,10:10:90)]
    "List of mean azimuth angles `[°]`"
    lazitab::Vector{FT} = collect(FT,5:10:355)

    # variables used for canopy_geometry!
    "Cosine of lazitab"
    cos_ttlo::Vector{FT} = cosd.(lazitab)
    "Cosine of lazitab - raa (relative azimuth angle), update with time"
    cos_philo::Vector{FT} = cosd.(lazitab .- 0)
    "Cosine of litab"
    cos_ttli::Vector{FT} = cosd.(litab)
    "Sine of litab"
    sin_ttli::Vector{FT} = sind.(litab)
    "Cache for volome scatter function"
    vol_scatt::Vector{FT} = ones(FT, 4)

    # This is changed afterwards, ignore here.
    "Inclination angles weight distribution"
    lidf::Vector{FT} = dladgen(LIDFa, LIDFb, litab_bnd)
    "List of level location (level = layer + 1)"
    xl::Vector{FT} = collect(FT, 0:-1.0/nLayer:-1)
    "1/nLayer"
    dx::FT = 1 / nLayer

    # local storage of dimension information
    "Number of azimuth angles"
    nAzi::Int = length(lazitab)
    "Number of inclination angles"
    nIncl::Int = length(litab)
end
