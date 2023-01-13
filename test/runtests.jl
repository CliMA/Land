using Land
using PkgUtility
using Test


pkgs = ["EmeraldConstants",
        "EarthSurface",
        "WaterPhysics",
        "EmeraldNamespace",
        "EmeraldOptics",
        "LeafOptics",
        "CanopyRadiativeTransfer",
        "Photosynthesis",
        "SoilHydraulics",
        "PlantHydraulics",
        "StomataModels",
        "SoilPlantAirContinuum"];

@testset verbose = true "CliMA Land Modules" begin
    for pkg in pkgs
        @info "Testing package $(pkg).jl...";
        include("../packages/$(pkg).jl/test/tests.jl");
    end;
end;
