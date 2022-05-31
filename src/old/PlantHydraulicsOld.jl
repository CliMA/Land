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

using ..PlantHydraulics: capacitance_buffer, critical_pressure, relative_hydraulic_conductance, xylem_end_pressure, xylem_pressure, xylem_pressure_profile!


# export public functions
export critical_flow, xylem_risk, plant_conductances!, xylem_flow


include("hydraulics/conductance.jl")
include("hydraulics/flow.jl"       )
include("hydraulics/capacitance.jl")


end # module
