using Test

@testset "Land" begin
  # include(joinpath(@__DIR__, "..", "experiments", "BRDF_test_directional.jl"))
  # include(joinpath(@__DIR__, "..", "experiments", "Canopy_netSWLW_APAR_fluxes.jl"))

  # include(joinpath(@__DIR__, "..", "experiments", "LeafPhoto_test.jl"))
  # include(joinpath(@__DIR__, "..", "experiments", "Optical_RTM_mSCOPE.jl"))
  # include(joinpath(@__DIR__, "..", "experiments", "RTM_source_test.jl"))
  # include(joinpath(@__DIR__, "..", "experiments", "run_field_data.jl"))
  # include(joinpath(@__DIR__, "..", "experiments", "step_light.jl"))
    include(joinpath(@__DIR__, "..", "experiments", "DiurnalCycle.jl"))
end

