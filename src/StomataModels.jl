module StomataModels

using ClimaCache: AbstractSoilVC, AirLayer, C3VJPModel, C4VJPModel, GCO₂Mode, Leaf
using ConstrainedRootSolvers: NewtonBisectionMethod, SolutionTolerance, find_zero
using DocStringExtensions: TYPEDFIELDS
using Photosynthesis: leaf_photosynthesis!, photosystem_temperature_dependence!
using PkgUtility: CP_D_MOL, GAS_R, K_STEFAN, M_H₂O, T_25
using PlantHydraulics: LeafHydraulics, MonoElementSPAC, critical_flow, relative_hydraulic_conductance, xylem_end_pressure
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack
using WaterPhysics: latent_heat_vapor, relative_diffusive_coefficient, relative_surface_tension, relative_viscosity, saturation_vapor_pressure


# export public types and structs
export BetaGLinearKleaf, BetaGLinearKsoil, BetaGLinearPleaf, BetaGLinearPsoil, BetaGLinearSWC, BetaVLinearKleaf, BetaVLinearKsoil, BetaVLinearPleaf, BetaVLinearPsoil, BetaVLinearSWC, CanopyLayer,
       EmpiricalStomatalModel, ESMBallBerry, ESMGentine, ESMLeuning, ESMMedlyn, GlcDrive, GswDrive, OptimizationStomatalModel, OSMEller, OSMSperry, OSMWang, OSMWAP, OSMWAPMod

# export public functions
export gsw_control!, gas_exchange!, prognostic_gsw!, solution_diff!, stomatal_conductance, update_leaf_AK!, update_leaf_TP!


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
include("model/risk.jl"       )
include("model/solution.jl"   )


end # module
