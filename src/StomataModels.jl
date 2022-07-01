module StomataModels

using ClimaCache: AbstractSoilVC, AbstractXylemVC, Leaf, Leaves1D, Leaves2D
using DocStringExtensions: METHODLIST
using PlantHydraulics: relative_hydraulic_conductance
using SoilHydraulics: relative_hydraulic_conductance
using UnPack: @unpack
using WaterPhysics: relative_diffusive_coefficient, relative_surface_tension


# include files
include("beta.jl"  )
include("limits.jl")


end # module
