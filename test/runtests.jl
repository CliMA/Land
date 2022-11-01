using Land
using Test


pkgs = ["EmeraldConstants", "WaterPhysics", "ClimaCache", "LeafOptics", "CanopyRadiativeTransfer", "Photosynthesis", "SoilHydraulics", "PlantHydraulics", "StomataModels", "SoilPlantAirContinuum"];

@testset verbose = true "CliMA Land Modules" begin
    for pkg in pkgs
        @info "Testing package $(pkg).jl...";
        include("../packages/$(pkg).jl/test/runtests.jl");
    end;
end;
