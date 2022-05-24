module PlantHydraulicsOld

using ClimaCache: AbstractPVCurve, AbstractSoilVC, AbstractXylemVC, BrooksCorey, ComplexVC, LinearPVCurve, LogisticVC, PowerVC, SegmentedPVCurve, VanGenuchten, WeibullVC
using ConstrainedRootSolvers: NewtonBisectionMethod, ReduceStepMethodND, SolutionTolerance, SolutionToleranceND, find_peak, find_zero
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using HypergeometricFunctions: _₂F₁
using PkgUtility: GAS_R, T_25, twarn, ρg_MPa
using SoilHydraulics: relative_hydraulic_conductance
using SpecialFunctions: gamma
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: relative_surface_tension, relative_viscosity

using ..PlantHydraulics: critical_pressure, relative_hydraulic_conductance, xylem_pressure

# export public types --- flow mode
export NonSteadyStateMode, SteadyStateMode

# export public types --- hydraulic tissue
export LeafHydraulics, RootHydraulics, StemHydraulics

# export public types --- hydraulic system
export GrassLikeOrganism, PalmLikeOrganism, TreeLikeOrganism, TreeSimple

# export public functions --- initialize plant
export create_grass, create_palm, create_tree

# export public functions
export flow_profile!, pressure_profile!, inititialize_legacy!, critical_flow, xylem_risk, plant_conductances!, roots_flow!, xylem_flow, update_PVF!, temperature_effects!, end_pressure


include("types/flow.jl" )
include("types/organ.jl")
include("types/plant.jl")

include("initialize/legacy.jl")
include("initialize/plant.jl" )

include("hydraulics/conductance.jl")
include("hydraulics/flow.jl"       )
include("hydraulics/pressure.jl"   )
include("hydraulics/temperature.jl")

include("hydraulics/capacitance.jl")


end # module
