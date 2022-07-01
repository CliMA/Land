module StomataModels

using ClimaCache: AbstractSoilVC, AbstractXylemVC
using DocStringExtensions: METHODLIST
using PlantHydraulics: relative_hydraulic_conductance
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack
using WaterPhysics: relative_surface_tension


# include files
include("beta.jl"  )
include("limits.jl")


end # module
