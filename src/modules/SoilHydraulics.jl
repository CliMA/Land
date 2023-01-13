module SoilHydraulics

using ConstrainedRootSolvers: ReduceStepMethodND, SolutionToleranceND, find_peak

using ..EmeraldConstants: CP_L, CP_L_MOL, M_H₂O, Λ_THERMAL_H₂O, ρ_H₂O, ρg_MPa
using ..WaterPhysics: relative_surface_tension, relative_viscosity
using ..EmeraldNamespace: MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC, NonSteadyStateFlow, Root, SteadyStateFlow, VanGenuchten

import ..EmeraldNamespace: BrooksCorey


include("../../packages/SoilHydraulics.jl/src/budget.jl"       )
include("../../packages/SoilHydraulics.jl/src/vulnerability.jl")


end # module
