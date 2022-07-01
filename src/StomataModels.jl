module StomataModels

using ClimaCache: AbstractSoilVC, BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ, LeafHydraulics, MonoElementSPAC, MonoMLGrassSPAC,
      MonoMLPalmSPAC, MonoMLTreeSPAC
using DocStringExtensions: METHODLIST
using PlantHydraulics: relative_hydraulic_conductance
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack
using WaterPhysics: relative_surface_tension


# include files
include("beta.jl")


end # module
