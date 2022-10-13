using CSV
using DataFrames
using Land
using Land.CanopyLayers
using Land.Photosynthesis
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using Pkg.Artifacts
using PkgUtility
using Test


ENV["JULIA_LOG_LEVEL"] = "WARN"


println("\nTesting the Land ODE functions...");
@testset "Land --- Land model" begin
    # initialize parameters
    @test true;
end
