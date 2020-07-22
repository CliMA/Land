module StomataModels

using CLIMAParameters
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters
using Photosynthesis
using PlantHydraulics
using WaterPhysics




# Define constants
Planet = CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH         = EarthParameterSet()
const K_25          = Planet.T_freeze(EARTH) + 25
const MOLMASS_WATER = Planet.molmass_water(EARTH)




# export public types and structs
export AbstractStomatalModel,
       CanopyLayer,
       EmpiricalStomatalModel,
       ESMBallBerry,
       ESMGentine,
       ESMLeuning,
       ESMMedlyn,
       OptimizationStomatalModel,
       OSMEller,
       OSMSperry,
       OSMWang,
       OSMWAP,
       OSMWAPMod




# export public functions
export empirical_gsw_from_model,
       envir_diff!,
       leaf_gsw_control!,
       leaf_photo_from_envir!,
       update_leaf_AK!,
       update_leaf_from_glc!,
       update_leaf_from_gsw!,
       update_leaf_TP!




include("types/canopylayer.jl"  )
include("types/stomatalmodel.jl")

include("empirical/general.jl"  )
include("empirical/ballberry.jl")
include("empirical/gentine.jl"  )
include("empirical/leuning.jl"  )
include("empirical/medlyn.jl"   )

include("optimal/general.jl")
include("optimal/eller.jl"  )
include("optimal/sperry.jl" )
include("optimal/wang.jl"   )
include("optimal/wap.jl"    )
include("optimal/wapmod.jl" )

include("gasexchange/control.jl")
include("gasexchange/refresh.jl")
include("gasexchange/update.jl" )




end # module
