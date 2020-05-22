module Land

using Parameters

include(joinpath("Utils", "PhysCon.jl"))
include(joinpath("Utils", "WaterVapor.jl"))

include(joinpath("Photosynthesis", "Photosynthesis.jl"))
include(joinpath("Radiation", "CanopyRT.jl"))

# module by Yujie, may need to merge with others
include("Photosynthesis/PhotosynthesisModels.jl")
include("Plant/Plant.jl"                        )

end # module
