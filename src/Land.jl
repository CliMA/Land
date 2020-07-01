module Land

# load constants and tools
include("Utils/LandParameters.jl"      )
include("Utils/Earth.jl"               )
include("Utils/MathTools.jl"           )
include("Utils/RootSolversExtension.jl")
include("Utils/WaterPhysics.jl"        )

# load the feature modules
include("Hydraulics/Hydraulics.jl"        )
include("Photosynthesis/Photosynthesis.jl")
include("Radiation/CanopyRT.jl"           )
include("SPAC/SPAC.jl"                    )

#include("Plant/Plant.jl")

end # module
