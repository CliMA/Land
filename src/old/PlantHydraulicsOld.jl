module PlantHydraulicsOld

using ClimaCache: AbstractPVCurve, AbstractSoilVC, AbstractXylemVC, BrooksCorey, ComplexVC, LeafHydraulics, LinearPVCurve, LogisticVC, MonoElementSPAC, MonoGrassSPAC, MonoPalmSPAC, MonoTreeSPAC,
      PowerVC, RootHydraulics, SegmentedPVCurve, StemHydraulics, VanGenuchten, WeibullVC
using ConstrainedRootSolvers: NewtonBisectionMethod, ReduceStepMethodND, SolutionTolerance, SolutionToleranceND, find_peak, find_zero
using DocStringExtensions: TYPEDEF, TYPEDFIELDS
using HypergeometricFunctions: _₂F₁
using PkgUtility: GAS_R, T_25, twarn, ρg_MPa
using SoilHydraulics: relative_hydraulic_conductance
using SpecialFunctions: gamma
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: relative_surface_tension, relative_viscosity

using ..PlantHydraulics: critical_pressure, relative_hydraulic_conductance, xylem_end_pressure, xylem_pressure


# export public functions
export flow_profile!, pressure_profile!, critical_flow, xylem_risk, plant_conductances!, roots_flow!, xylem_flow, update_PVF!


include("hydraulics/conductance.jl")
include("hydraulics/flow.jl"       )
include("hydraulics/pressure.jl"   )
include("hydraulics/capacitance.jl")


end # module
