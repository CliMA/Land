using Test

ENV["JULIA_LOG_LEVEL"] = "WARN"

@testset "Land" begin
    include(joinpath(@__DIR__, "..", "experiments", "DiurnalCycle.jl"))
    include(joinpath(@__DIR__, "..", "experiments", "Radiation_Test_BRDF.jl"))
end

# test the plant module isnan and FT consistency
include("test_plant.jl")