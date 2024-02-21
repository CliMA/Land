###############################################################################
#
# SoilPlantAirContinuum system for mono species system
#
###############################################################################
"""
    mutable struct SPACMono{FT}

Struct that mono species SoilPlantAirContinuum system.

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct SPACMono{FT<:AbstractFloat}
    "Input file for SPAC spectra used in wl_set and soil_opt"
    opti_file::String = LAND_2021

    "Soil layers bounds `[m]`"
    soil_bounds::Vector{FT} = FT[0,-0.1,-0.2,-0.3,-0.5,-0.8,-1.2,-2.0]
    "Air layers bounds `[m]`"
    air_bounds ::Vector{FT} = collect(FT,0:1:20)

    "Root depth `[m]`"
    z_root  ::FT = FT(-1)
    "Canopy maximal height `[m]`"
    z_canopy::FT = FT(10)

    "Plant hydraulic system"
    plant_hs::Union{GrassLikeOrganism{FT}, PalmLikeOrganism{FT}, TreeLikeOrganism{FT}} = create_grass(z_root, z_canopy,soil_bounds, air_bounds)
    "Number of canopy layers"
    n_canopy::Int = length(plant_hs.canopy_index_in_air)
    "Number of root layers"
    n_root  ::Int = length(plant_hs.root_index_in_soil)
    "Plant photosynthesis systems"
    plant_ps::Vector{CanopyLayer{FT}} = [CanopyLayer{FT}() for i in 1:n_canopy]
    "Basal area `[m²]`"
    ba::FT = sum( [plant_hs.roots[i].area for i in 1:n_root] )
    "Ground area `[m²]`"
    ga::FT = ba * 500
    "Leaf area `[m²]`"
    la::FT = sum( plant_hs.leaves[i].area for i in 1:n_canopy )

    # Surrounding AirLayer
    "Air layers"
    envirs::Vector{AirLayer{FT}} = [AirLayer{FT}() for i in 1:n_canopy]

    # Wind related
    "Aerodynamic roughness `[m]`"
    wind_z0::FT = z_canopy * FT(0.07)
    "Zero plane displacement `[m]`"
    wind_d ::FT = z_canopy * 2/3
    "Mean layer height `[m]`"
    wind_zs::Vector{FT} = [(air_bounds[i]+air_bounds[i+1])/2 for i in 1:length(air_bounds)-1]
    "Wind speed per layer `[m s⁻¹]`"
    winds::Vector{FT} = [FT(1) for i in 1:length(air_bounds)-1]

    # Soil layers information
    # TODO bridge Soil module later
    "Maximal soil water content"
    mswc  ::Vector{FT} = [plant_hs.roots[i].sh.Θs for i in 1:n_root]
    "Current soil water content"
    swc   ::Vector{FT} = [plant_hs.roots[i].sh.Θs for i in 1:n_root]
    "Array of soil matric potential `[MPa]`"
    p_soil::Vector{FT} = FT[0 for i in 1:n_root]
    "Maximal soil depth `[m]`"
    h_soil::FT = abs(z_root)

    # geography related
    "Latitude `[°]`"
    latitude ::FT = 30
    "Longitude `[°]`"
    longitude::FT = 116
    "Elevation `[m]`"
    elevation::FT = 0

    # photosynthesis mode and stomatal model scheme
    "Photosynthesis parameter set"
    photo_set::Union{C3ParaSet{FT}, C4ParaSet{FT}} = C3CLM(FT)
    "Stomatal behavior scheme"
    stomata_model::Union{EmpiricalStomatalModel{FT}, OptimizationStomatalModel{FT}} = ESMBallBerry{FT}()

    # For CanopyLayers module
    "Solar angle container"
    angles::SolarAngles{FT} = SolarAngles{FT}()
    "Canopy4RT container"
    canopy_rt::Canopy4RT{FT} = Canopy4RT{FT}(nLayer = n_canopy, LAI = la/ga)
    "Wave length container"
    wl_set::WaveLengths{FT} = WaveLengths{FT}(opti_file = opti_file)
    "RT dimensions"
    rt_dim::RTDimensions = RTDimensions(canopy_rt, wl_set)
    "CanopyRads container"
    can_rad::CanopyRads{FT} = CanopyRads{FT}(rt_dim)
    "CanopyOpticals container"
    can_opt::CanopyOpticals{FT} = CanopyOpticals{FT}(rt_dim)
    "Array of LeafBios container"
    leaves_rt::Vector{LeafBios{FT}} = [LeafBios{FT}(rt_dim) for i in 1:n_canopy]
    "SoilOpticals container"
    soil_opt::SoilOpticals{FT} = SoilOpticals(wl_set)
    "Incoming radiation container"
    in_rad::IncomingRadiation{FT} = IncomingRadiation(wl_set)
    "RT container"
    rt_con::RTCache{FT} = RTCache{FT}(rt_dim)
    "Container for sunlit leaf area fraction in each layer"
    f_SL::Vector{FT} = repeat(canopy_rt.lidf, outer=[ canopy_rt.nAzi ]) / canopy_rt.nAzi;

    # local storage for canopy GPP and NPP
    "Canopy GPP per ground area"
    f_gpp::FT = 0
    "Canopy GPP per ground area"
    f_npp::FT = 0
    "Canopy water flux per ground area"
    f_H₂O::FT = 0
end
