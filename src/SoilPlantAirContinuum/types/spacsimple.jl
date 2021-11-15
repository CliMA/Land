###############################################################################
#
# SoilPlantAirContinuum system for simple mono species system
#
###############################################################################
"""
    mutable struct SPACSimple{FT}

Struct that simplest mono species SoilPlantAirContinuum system, with 1 root,
    stem, and leaf.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SPACSimple{FT<:AbstractFloat}
    # Hydraulic system
    "Hydraulic systems"
    hs::TreeSimple{FT} = TreeSimple{FT}()
    "Critical flow rate"
    ec::FT = FT(50);

    # Photosynthesis parameters
    "Photosynthesis system"
    ps::Leaf{FT} = Leaf{FT}()
    "Ratio between Vcmax25 and Jmax25"
    vtoj::FT = FT(1.67)

    # Surrounding AirLayer
    "Environmental conditions"
    envir::AirLayer{FT} = AirLayer{FT}()

    # local container for returned results
    "Container for gas exchange for a layer"
    container1L::SPACContainer1L{FT} = SPACContainer1L{FT}()
    "Container for gas exchange of sunlit and shaded layers"
    container2L::SPACContainer2L{FT} = SPACContainer2L{FT}()
    "Container for default hydraulic conductance"
    containerKS::Array{FT,1} = FT[hs.root.k_max, hs.stem.k_max, hs.leaf.k_sla]
    "Container for optimizer"
    containerOP::FT = FT(0)
    "Container for optimal sunlit layer flow rate"
    opt_f_sl   ::FT = FT(0)
    "Container for optimal shaded layer flow rate"
    opt_f_sh   ::FT = FT(0)
    "Container for optimal leaf area per basal area"
    opt_laba   ::FT = FT(2000)
    "Container for optimal Vcmax25"
    opt_vmax   ::FT = FT(5)

    # leaf related
    "Leaf area index"
    lai  ::FT = FT(3)
    "Leaf area per basal area"
    laba ::FT = 1500
    "Maximal stomatal conductance limit at 25 °C"
    g_max::FT = 0.8
    "Ground area per basal area"
    gaba ::FT = 500
    "Leaf width"
    width::FT = 0.05

    # soil related
    "Maximal soil water content"
    mswc  ::FT = (hs.root.sh).Θs
    "Current soil water content"
    swc   ::FT = (hs.root.sh).Θs
    "Soil matrical water potential"
    p_soil::FT = 0.0
    "Soil depth, 2X mean root depth"
    h_soil::FT = 2

    # leaf invest related
    "Leaf construction cost per leaf area"
    c_cons::FT = 2
    "Leaf nutrient cost per Vcmax25 per leaf area"
    c_vmax::FT = 0.04

    # geography related
    "Latitude `[°]`"
    latitude ::FT = 30
    "Longitude `[°]`"
    longitude::FT = 116
    "Elevation `[m]`"
    elevation::FT = 0
end
