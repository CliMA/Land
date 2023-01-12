module PlantHydraulics

using ConstrainedRootSolvers: NewtonBisectionMethod, SolutionTolerance, find_zero
using Statistics: mean
using UnPack: @unpack

#using EmeraldConstants: CP_D_MOL, CP_L_MOL, GAS_R, M_H₂O, T₂₅, ρg_MPa
#using WaterPhysics: latent_heat_vapor, relative_surface_tension, relative_viscosity, saturation_vapor_pressure
#using EmeraldNamespace: AbstractSoilVC, AbstractXylemVC, ComplexVC, LinearPVCurve, LogisticVC, PowerVC, SegmentedPVCurve, WeibullVC
#using EmeraldNamespace: BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ
#using EmeraldNamespace: AbstractStomataModel, AndereggSM, BallBerrySM, EllerSM, GentineSM, LeuningSM, MedlynSM, SperrySM, WangSM, Wang2SM
#using EmeraldNamespace: Leaf, LeafHydraulics, Leaves1D, Leaves2D, NonSteadyStateFlow, Root, RootHydraulics, Soil, SoilLayer, SteadyStateFlow, Stem, StemHydraulics
#using EmeraldNamespace: MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
#using SoilHydraulics: soil_θ, soil_ψ_25

#import SoilHydraulics: relative_hydraulic_conductance


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
