module PlantHydraulics

using ClimaCache: ComplexVC, Leaf, LeafHydraulics, Leaves1D, Leaves2D, LinearPVCurve, LogisticVC, MonoElementSPAC, MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC, NonSteadyStateFlow, PowerVC, Root,
      RootHydraulics, SegmentedPVCurve, SteadyStateFlow, Stem, StemHydraulics, WeibullVC
using ClimaCache: GAS_R, T_25, ρg_MPa
using ConstrainedRootSolvers: NewtonBisectionMethod, SolutionTolerance, find_zero
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: relative_surface_tension, relative_viscosity

import SoilHydraulics: relative_hydraulic_conductance


# include functions
include("critical_pressure.jl")
include("derivative.jl"       )
include("flow_profile.jl"     )
include("legacy.jl"           )
include("pressure_profile.jl" )
include("pressure_volume.jl"  )
include("target_flow.jl"      )
include("vulnerability.jl"    )


end # module
