module Land

# load constants and tools
include("Utils/LandParameters.jl")
include("Utils/MathTools.jl")
include("Utils/WaterPhysics.jl")

# load the feature modules
include("Photosynthesis/Photosynthesis.jl")

include("Plant/Plant.jl")
include("Leaf/Leaf.jl")
include("Radiation/CanopyRT.jl")
include("SPAC/SPAC.jl")

end # module
