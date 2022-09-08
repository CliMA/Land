using ClimaCache
using SoilPlantAirContinuum
using Test


@testset verbose = true "SoilPlantAirContinuum Test" begin
    @testset "Core Function" begin
        for FT in [Float32, Float64]
            for spac in [ClimaCache.MonoMLGrassSPAC{FT}(), ClimaCache.MonoMLPalmSPAC{FT}(), ClimaCache.MonoMLTreeSPAC{FT}()]
                SoilPlantAirContinuum.initialize!(spac);
                @test true;
                SoilPlantAirContinuum.soil_plant_air_continuum!(spac, FT(1));
                @test true;
                SoilPlantAirContinuum.soil_plant_air_continuum!(spac);
                @test true;
            end;
        end;
    end;

    @testset "Quantities" begin
        for FT in [Float32, Float64]
            for spac in [ClimaCache.MonoMLGrassSPAC{FT}(), ClimaCache.MonoMLPalmSPAC{FT}(), ClimaCache.MonoMLTreeSPAC{FT}()]
                SoilPlantAirContinuum.initialize!(spac);
                SoilPlantAirContinuum.soil_plant_air_continuum!(spac, FT(1));
                @test !isnan(SoilPlantAirContinuum.gross_primary_productivity(spac));
                @test !isnan(SoilPlantAirContinuum.transpiration_rate(spac));
            end;
        end;
    end;
end;
