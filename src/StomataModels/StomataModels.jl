module StomataModels

using CLIMAParameters
using CLIMAParameters.Planet
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters
using ..Photosynthesis
using ..PlantHydraulics
using WaterPhysics




# Define constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH        = EarthParameterSet();
CP_DRYAIR(FT)      = FT( cp_d(EARTH) );
K_0(FT)            = FT( T_freeze(EARTH) );
K_STEFAN(FT)       = FT( Stefan() );
MOLMASS_DRYAIR(FT) = FT( molmass_dryair(EARTH) );
MOLMASS_WATER(FT)  = FT( molmass_water(EARTH) );

K_25(FT)           = K_0(FT) + 25;
CP_DRYAIR_MOL(FT)  = CP_DRYAIR(FT) * MOLMASS_DRYAIR(FT);




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
export gsw_control!,
       gas_exchange!,
       prognostic_gsw!,
       solution_diff!,
       stomatal_conductance,
       update_leaf_AK!,
       update_leaf_TP!




include("types/beta.jl"         )
include("types/canopylayer.jl"  )
include("types/drive.jl"        )
include("types/stomatalmodel.jl")

include("model/beta.jl"       )
include("model/control.jl"    )
include("model/empirical.jl"  )
include("model/gasexchange.jl")
include("model/nocturnal.jl"  )
include("model/prognostic.jl" )
include("model/refresh.jl"    )
include("model/solution.jl"   )




end # module
