module PlantHydraulics

using CLIMAParameters
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters




# define global constants
Planet = CLIMAParameters.Planet
struct EarthParameterSet <: AbstractEarthParameterSet end
const EARTH   = EarthParameterSet()
const GRAVITY = Planet.grav(EARTH)
const K_25    = Planet.T_freeze(EARTH) + 25
const ρ_H₂O   = Planet.ρ_cloud_liq(EARTH)
const ρg_MPa  = ρ_H₂O * GRAVITY * 1e-6




# export public types
export AbstractHydraulicSystem,
       LeafHydraulics,
       RootHydraulics,
       StemHydraulics,
       AbstractPlantHS,
       GrassLikeHS,
       PalmLikeHS,
       TreeLikeHS,
       AbstractVulnerability,
       WeibullDual,
       WeibullSingle

# export public functions
export create_grass_like_hs,
       create_palm_like_hs,
       create_tree_like_hs,
       hydraulic_p_profile!,
       leaf_e_crit,
       leaf_xylem_risk,
       q_diff,
       root_qs_p_from_q,
       tree_e_crit,
       tree_p_from_flow,
       xylem_k_ratio,
       xylem_p_from_flow




include("types/vulnerability.jl"   )
include("types/hydraulics_organ.jl")
include("types/hydraulics_plant.jl")
include("types/initialize_plant.jl")
include("vulnerability/weibull.jl" )
include("hydraulics/base.jl"       )
include("hydraulics/leaf.jl"       )
include("hydraulics/plant.jl"      )
include("hydraulics/root.jl"       )




end # module
