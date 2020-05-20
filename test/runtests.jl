using Test

ENV["JULIA_LOG_LEVEL"] = "WARN"


@testset "Land" begin
    include(joinpath(@__DIR__, "..", "experiments", "DiurnalCycle.jl"))
    include(joinpath(@__DIR__, "..", "experiments", "Radiation_Test_BRDF.jl"))
end

# test the plant module isnan and FT consistency
include("Plants/test_nan_and_FT.jl")

# test the plant-rt interface
include("Plants/test_rt_interface.jl")
