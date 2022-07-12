module SoilPlantAirContinuumNew

using ClimaCache: AirLayer, MonoMLTreeSPAC
using WaterPhysics: saturation_vapor_pressure


include("model.jl" )
include("update.jl")


end # module
