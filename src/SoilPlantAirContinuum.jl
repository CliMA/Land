module SoilPlantAirContinuum

using ClimaCache: AirLayer


include("environment.jl")

# the archive functions from older version
include("SoilPlantAirContinuumOld.jl")

end # module
