module PlantHydraulics

using ClimaCache: AbstractSoilVC, AbstractXylemVC, AndereggSM, BallBerrySM, BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ, ComplexVC,
      EllerSM, GentineSM, Leaf, LeafHydraulics, Leaves1D, Leaves2D, LeuningSM, LinearPVCurve, LogisticVC, MedlynSM, MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC,
      NonSteadyStateFlow, PowerVC, Root, RootHydraulics, SegmentedPVCurve, SperrySM, SteadyStateFlow, Stem, StemHydraulics, WangSM, Wang2SM, WeibullVC
using ClimaCache: GAS_R, T_25, ρg_MPa
using ConstrainedRootSolvers: NewtonBisectionMethod, SolutionTolerance, find_zero
using SoilHydraulics: soil_θ, soil_ψ_25
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: relative_surface_tension, relative_viscosity

import SoilHydraulics: relative_hydraulic_conductance


# include functions
include("beta.jl"             )
include("critical_pressure.jl")
include("derivative.jl"       )
include("flow_profile.jl"     )
include("legacy.jl"           )
include("pressure_profile.jl" )
include("pressure_volume.jl"  )
include("target_flow.jl"      )
include("vulnerability.jl"    )


end # module
