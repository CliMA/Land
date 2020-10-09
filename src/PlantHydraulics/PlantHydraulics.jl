module PlantHydraulics

using CLIMAParameters
using CLIMAParameters.Planet
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters
using Statistics
using WaterPhysics




# define global constants
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH = EarthParameterSet()
GRAVITY(FT) = FT( grav(EARTH) );
K_25(FT)    = FT( T_freeze(EARTH) ) + 25;
ρ_H₂O(FT)   = FT( ρ_cloud_liq(EARTH) );
ρg_MPa(FT)  = ρ_H₂O(FT) * GRAVITY(FT) * FT(1e-6);




# export public types
export AbstractXylemVC,
       WeibullDual,
       WeibullSingle,
       AbstractSoilVC,
       BrooksCorey,
       VanGenuchten,
       AbstractHydraulicSystem,
       LeafHydraulics,
       RootHydraulics,
       StemHydraulics,
       AbstractPlantHS,
       GrassLikeHS,
       PalmLikeHS,
       TreeLikeHS,
       TreeSimple

# export public functions
export create_grass_like_hs,
       create_palm_like_hs,
       create_tree_like_hs,
       hydraulic_p_profile!,
       inititialize_legacy!,
       leaf_e_crit,
       leaf_xylem_risk,
       plant_conductances!,
       recalculate_roots_flow!,
       root_q_from_pressure,
       roots_flow!,
       soil_erwc,
       soil_k_ratio_erwc,
       soil_k_ratio_p25,
       soil_k_ratio_rwc,
       soil_k_ratio_swc,
       soil_p_25_erwc,
       soil_p_25_rwc,
       soil_p_25_swc,
       soil_rwc,
       soil_swc,
       tree_e_crit,
       vc_temperature_effects!,
       xylem_k_ratio,
       xylem_p_crit,
       xylem_p_from_flow




include("types/vulnerability.jl"   )
include("types/hydraulics_organ.jl")
include("types/hydraulics_plant.jl")
include("types/initialize_plant.jl")

include("vulnerability/pcrit.jl")
include("vulnerability/soil.jl" )
include("vulnerability/vc.jl"   )

include("hydraulics/leaf.jl"       )
include("hydraulics/legacy.jl"     )
include("hydraulics/plant.jl"      )
include("hydraulics/pressure.jl"   )
include("hydraulics/root.jl"       )
include("hydraulics/temperature.jl")




end # module
