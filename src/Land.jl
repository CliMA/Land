module Land

using Revise


# include sub modules
include("CanopyLayers/CanopyLayers.jl")
include("Photosynthesis/Photosynthesis.jl")
include("PlantHydraulics/PlantHydraulics.jl")
include("StomataModels/StomataModels.jl")
include("SoilPlantAirContinuum/SoilPlantAirContinuum.jl")


end # module
