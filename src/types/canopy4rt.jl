###############################################################################
#
# Canopy for RT calculations
#
###############################################################################
"""
    struct Canopy4RT

A canopy struct for the radiation transfer module

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct Canopy4RT{FT<:AbstractFloat}
    "Number of canopy layers"
    nLayer::Int = 20
    "Leaf Area Index"
    LAI       ::FT  = FT(3.0 )
    "Clumping factor"
    Ω         ::FT  = FT(1.0 )
    "Structure factor a"
    clump_a   ::FT  = FT(1.0 )
    "Structure factor b"
    clump_b   ::FT  = FT(0.0 )
    "Leaf width"
    leaf_width::FT  = FT(0.1 )
    "Vegetation height"
    hc        ::FT  = FT(2.0 )
    "Leaf Inclination"
    LIDFa     ::FT  = FT(0.0 )
    "Variation in leaf inclination"
    LIDFb     ::FT  = FT(0.0 )
    "HotSpot parameter (still need to check!)"
    hot       ::FT  = FT(0.05)

    # tree/canopy/leaf traits
    "Canopy height `[m]`"
    height::FT = FT(20.0  )
    "Canopy roughness `[m]`"
    z0m   ::FT = FT(1.0   )
    "Tree roughtnes `[m]`"
    z0h   ::FT = FT(-999.0)
    "Canopy displacement height `[m]`"
    d     ::FT = FT(-999.0)
    "m/sqrt(s) turbulent transfer coefficient"
    Cd    ::FT = FT(0.01  )


    # Some more derived parameters:
    "List of mean inclination angles `[°]`"
    litab    ::Array{FT,1} = collect(FT,5:10:85)
    "List of inclination angle boundaries `[°]`"
    litab_bnd::Array{FT,2} = [collect(0:10:80) collect(FT,10:10:90)]
    "List of mean azimuth angles `[°]`"
    lazitab  ::Array{FT,1} = collect(FT,5:10:355)

    # This is changed afterwards, ignore here.
    "Inclination angles weight distribution"
    lidf::Array{FT,1} = litab .* 0
    "List of level location (level = layer + 1)"
    xl  ::Array{FT,1} = collect(0.0:-1.0/nLayer:-1.0)
    "1/nlayers"
    dx  ::FT          = FT(1)/nLayer
end
