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




# TODO add extra-xylary VC
# TODO add pressure-volume curve functions




# export public types
export AbstractHydraulicSystem,
       LeafHydraulics,
       RootHydraulics,
       Palm,
       StemHydraulics,
       Tree

# export public functions
export create_palm,
       create_tree,
       hydraulic_p_profile!,
       leaf_e_crit,
       leaf_xylem_risk,
       root_qs_p_from_q,
       tree_e_crit,
       tree_p_from_flow,
       weibull_k_ratio,
       xylem_p_from_flow




include("types.jl"  )
include("testing.jl")
include("weibull.jl")
include("flow.jl"   )




end # module
