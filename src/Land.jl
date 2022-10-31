module Land


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
