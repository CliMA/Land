module StomataModels

using CLIMAParameters
using CLIMAParameters.Planet
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters
using Photosynthesis
using PlantHydraulics
using WaterPhysics




# Define constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH       = EarthParameterSet()
K_0(FT)           = FT( T_freeze(EARTH) );
K_25(FT)          = K_0(FT) + 25;
MOLMASS_WATER(FT) = FT( molmass_water(EARTH) );




# export public types and structs
export AbstractBetaFunction,
       AbstractBetaG,
       AbstractBetaV,
       AbstractDrive,
       AbstractStomatalModel,
       BetaGLinearKleaf,
       BetaGLinearPleaf,
       BetaGLinearPsoil,
       BetaGLinearSWC,
       BetaVLinearKleaf,
       BetaVLinearPleaf,
       BetaVLinearPsoil,
       BetaVLinearSWC,
       CanopyLayer,
       EmpiricalStomatalModel,
       ESMBallBerry,
       ESMGentine,
       ESMLeuning,
       ESMMedlyn,
       GlcDrive,
       GswDrive,
       OptimizationStomatalModel,
       OSMEller,
       OSMSperry,
       OSMWang,
       OSMWAP,
       OSMWAPMod

# export public functions
export leaf_gsw_control!,
       gas_exchange!,
       solution_diff!,
       stomatal_conductance,
       update_leaf_AK!,
       update_leaf_TP!




include("types/beta.jl"         )
include("types/canopylayer.jl"  )
include("types/drive.jl"        )
include("types/stomatalmodel.jl")

include("empirical/beta.jl"     )
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
