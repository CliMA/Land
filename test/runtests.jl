using Land
using Pkg
using Test


@testset verbose = true "CliMA Land" begin
    for pkgname in [
                "EmeraldConstants",
                "WaterPhysics",
                "ClimaCache",
                "LeafOptics",
                "CanopyRadiativeTransfer",
                "Photosynthesis",
                "SoilHydraulics",
                "PlantHydraulics",
                "StomataModels",
                "SoilPlantAirContinuum"]
        Pkg.test(pkgname);
        @test true;
    end;
end;
