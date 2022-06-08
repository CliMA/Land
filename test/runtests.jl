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
            CanopyRadiativeTransfer.canopy_geometry!(can, angles);
            @test true;
        end;
    end;
end;
