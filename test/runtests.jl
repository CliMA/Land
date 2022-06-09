using CanopyRadiativeTransfer
using ClimaCache
using Test


@testset verbose = true "CanopyRadiativeTransfer Test" begin
    @testset "Leaf Inclination Angles" begin
        for FT in [Float32, Float64]
            can = ClimaCache.HyperspectralMLCanopy{FT}();
            CanopyRadiativeTransfer.inclination_angles!(can, can.LIDF, FT(0), FT(0));
            @test true;
        end;
    end;

    @testset "Clumping Index" begin
        for FT in [Float32, Float64]
            can = ClimaCache.HyperspectralMLCanopy{FT}();
            angles = ClimaCache.SunSensorGeometry{FT}();
            CanopyRadiativeTransfer.clumping_index!(can, angles);
            @test true;
        end;
    end;

    @testset "Hyperspectral RT" begin
        for FT in [Float32, Float64]
            can = ClimaCache.HyperspectralMLCanopy{FT}();
            angles = ClimaCache.SunSensorGeometry{FT}();
            leaf = ClimaCache.Leaf{FT}("C3");
            leaves = [deepcopy(leaf) for i in 1:20];
            soil = ClimaCache.Soil{FT}(FT[0,-1]);
            rad = ClimaCache.HyperspectralRadiation{FT}();
            CanopyRadiativeTransfer.canopy_optical_properties!(can, angles);
            @test true;
            CanopyRadiativeTransfer.canopy_optical_properties!(can, leaves, soil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(can, rad);
            @test true;
        end;
    end;
end;
