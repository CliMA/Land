module PlantHydraulics

using CLIMAParameters
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters
using WaterPhysics




# define global constants
Planet = CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH   = EarthParameterSet()
const GRAVITY = Planet.grav(EARTH)
const K_25    = Planet.T_freeze(EARTH) + 25
const ρ_H₂O   = Planet.ρ_cloud_liq(EARTH)
const ρg_MPa  = ρ_H₂O * GRAVITY * 1e-6




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
       q_diff,
       root_qs_p_from_q,
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
