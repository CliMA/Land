module Land

using DiffEqOperators
using DocStringExtensions
using OrdinaryDiffEq
using Parameters
using PkgUtility: T_25
using WaterPhysics: TraceGasAir, TraceGasCO₂, diffusive_coefficient

# export public types
export VerticalLayers

# export public functions
export vertical_layers!

# include sub modules
include("CanopyLayers/CanopyLayers.jl"                  )
include("Photosynthesis/Photosynthesis.jl"              )
include("PlantHydraulics/PlantHydraulics.jl"            )
include("StomataModels/StomataModels.jl"                )
include("SoilPlantAirContinuum/SoilPlantAirContinuum.jl")

# include types
include("Land/Types.jl"    )
include("Land/AirLayers.jl")

end # module
