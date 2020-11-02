using CSV
using DataFrames
using Land
using Land.CanopyLayers
using Land.Photosynthesis
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using Pkg.Artifacts
using Test

ENV["JULIA_LOG_LEVEL"] = "WARN"




include("recursive.jl")




include("test_CanopyLayers.jl"   )
include("test_Photosynthesis.jl" )
include("test_PlantHydraulics.jl")
include("test_StomataModels.jl"  )
include("test_SPAC.jl"           )
include("test_Land.jl"           )
