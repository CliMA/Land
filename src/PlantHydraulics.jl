module PlantHydraulics

using CLIMAParameters
using ConstrainedRootSolvers
using DocStringExtensions
using Parameters




# define global constants
Planet = CLIMAParameters.Planet
""" Struct inherited from AbstractEarthParameterSet """
struct EarthParameterSet <: AbstractEarthParameterSet end
""" Earth as EarthParameterSet """
const EARTH   = EarthParameterSet()
""" Gravitational acceleration `[m s⁻²]` """
const GRAVITY = Planet.grav(EARTH)
""" Temperature at 25 Celcius `[K]` """
const K_25    = Planet.T_freeze(EARTH) + 25
""" Density of water `[kg m⁻³]`"""
const ρ_H₂O   = Planet.ρ_cloud_liq(EARTH)
""" Product of ρ and g `[MPa m⁻¹]` """
const ρg_MPa  = ρ_H₂O * GRAVITY * 1e-6




# TODO add extra-xylary VC
# TODO add pressure-volume curve functions




# export public types
export AbstractHydraulicSystem,
       LeafHydraulics,
       RootHydraulics,
       StemHydraulics,
       Tree

# export public functions
export create_tree,
       hydraulic_p_profile!,
       leaf_e_crit,
       leaf_xylem_risk,
       root_qs_p_from_q,
       weibull_k_ratio,
       xylem_p_from_flow




include("types.jl"  )
include("testing.jl")
include("weibull.jl")
include("flow.jl"   )




end # module
