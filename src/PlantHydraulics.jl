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




# export public types --- soil vulnerability
export AbstractSoilVC,
       BrooksCorey,
       VanGenuchten

# export public types --- xylem vulnerability
export AbstractXylemVC,
       WeibullDual,
       WeibullSingle

# export public types --- pressure volume curve
export AbstractCapacity,
       PVCurveLinear

# export public types --- flow mode
export AbstractFlowMode,
       NonSteadyStateMode,
       SteadyStateMode

# export public types --- hydraulic tissue
export AbstractHydraulicSystem,
       LeafHydraulics,
       RootHydraulics,
       StemHydraulics

# export public types --- hydraulic system
export AbstractPlantHS,
       GrassLikeHS,
       PalmLikeHS,
       TreeLikeHS,
       TreeSimple

# export public functions --- initialize soil type
export create_VanGenuchten,
       VGClay,
       VGClayLoam,
       VGLoamySand,
       VGLoamySand2,
       VGSand,
       VGSandyClayLoam,
       VGSandyClayLoam2,
       VGSandyLoam,
       VGSilt,
       VGSiltyClay,
       VGSiltyClayLoam,
       VGSiltyLoam

# export public functions --- initialize plant
export create_grass_like_hs,
       create_palm_like_hs,
       create_tree_like_hs

# export public functions --- curves related
export p_from_volume,
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
       xylem_k_ratio,
       xylem_p_crit

# export public functions
export hydraulic_p_profile!,
       inititialize_legacy!,
       leaf_e_crit,
       leaf_xylem_risk,
       plant_conductances!,
       recalculate_roots_flow!,
       root_q_from_pressure,
       roots_flow!,
       tree_e_crit,
       update_PVF!,
       vc_temperature_effects!,
       xylem_pressure




include("types/curves.jl")
include("types/flow.jl"  )
include("types/tissue.jl")
include("types/plant.jl" )

include("initialize/soil.jl" )
include("initialize/plant.jl")

include("curves/capacity.jl")
include("curves/soil.jl"    )
include("curves/xylem.jl"   )

include("hydraulics/xylem.jl")

include("hydraulics/capacitance.jl")
include("hydraulics/leaf.jl"       )
include("hydraulics/legacy.jl"     )
include("hydraulics/plant.jl"      )
include("hydraulics/pressure.jl"   )
include("hydraulics/root.jl"       )
include("hydraulics/temperature.jl")




end # module
