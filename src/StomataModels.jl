module StomataModels

using ClimaCache: AbstractSoilVC, BetaFunction, BetaParameterKleaf, BetaParameterKsoil, BetaParameterPleaf, BetaParameterPsoil, BetaParameterΘ, LeafHydraulics
using DocStringExtensions: METHODLIST
using PlantHydraulics: relative_hydraulic_conductance
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack


# include files
include("beta.jl")


end # module
