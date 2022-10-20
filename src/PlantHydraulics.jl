module PlantHydraulics

using ClimaCache: AbstractSoilVC, AbstractXylemVC, AndereggSM, BallBerrySM, BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ, ComplexVC,
      EllerSM, GentineSM, Leaf, LeafHydraulics, Leaves1D, Leaves2D, LeuningSM, LinearPVCurve, LogisticVC, MedlynSM, MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC,
      NonSteadyStateFlow, PowerVC, Root, RootHydraulics, SegmentedPVCurve, Soil, SoilLayer, SperrySM, SteadyStateFlow, Stem, StemHydraulics, WangSM, Wang2SM, WeibullVC
using ConstrainedRootSolvers: NewtonBisectionMethod, SolutionTolerance, find_zero
using EmeraldConstants: CP_D_MOL, CP_L_MOL, GAS_R, M_H₂O, T₂₅, ρg_MPa
using SoilHydraulics: soil_θ, soil_ψ_25
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure

import SoilHydraulics: relative_hydraulic_conductance


# include functions
include("beta.jl"             )
include("budget.jl"           )
include("critical_pressure.jl")
include("derivative.jl"       )
include("flow_profile.jl"     )
include("legacy.jl"           )
include("pressure_profile.jl" )
include("pressure_volume.jl"  )
include("target_flow.jl"      )
include("vulnerability.jl"    )


end # module
