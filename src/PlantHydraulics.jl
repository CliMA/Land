module PlantHydraulics

using ClimaCache: ComplexVC, Leaf, LeafHydraulics, LinearPVCurve, LogisticVC, MonoElementSPAC, MonoGrassSPAC, MonoPalmSPAC, MonoTreeSPAC, NonSteadyStateFlow, PowerVC, Root, RootHydraulics,
      SegmentedPVCurve, SteadyStateFlow, Stem, StemHydraulics, WeibullVC
using DocStringExtensions: METHODLIST
using PkgUtility: GAS_R, T_25, ρg_MPa
using Statistics: mean
using UnPack: @unpack
using WaterPhysics: relative_surface_tension, relative_viscosity

import SoilHydraulics: relative_hydraulic_conductance


# export public types from ClimaCache
export ComplexVC, LeafHydraulics, LinearPVCurve, LogisticVC, MonoElementSPAC, MonoGrassSPAC, MonoPalmSPAC, MonoTreeSPAC, PowerVC, RootHydraulics, SegmentedPVCurve, StemHydraulics, WeibullVC

# export public functions
export xylem_flow_profile!, xylem_pressure_profile!


# include functions
include("critical_pressure.jl")
include("flow_profile.jl"     )
include("legacy.jl"           )
include("pressure_profile.jl" )
include("pressure_volume.jl"  )
include("vulnerability.jl"    )


end # module
