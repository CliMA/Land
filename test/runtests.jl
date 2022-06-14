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
            hsoil = ClimaCache.Soil{FT}(FT[0,-1]);
            bsoil = ClimaCache.Soil{FT}(FT[0,-1]; n_λ = 1);
            rad = ClimaCache.HyperspectralRadiation{FT}();
            CanopyRadiativeTransfer.canopy_optical_properties!(can, angles);
            @test true;
            CanopyRadiativeTransfer.canopy_optical_properties!(can, leaves, hsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_optical_properties!(can, leaves, bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(can, leaves, rad, hsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(can, leaves, rad, bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(can, leaves, FT(100), hsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(can, leaves, FT(100), bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_fluorescence!(can, leaves);
            @test true;
        end;
    end;

    @testset "Remote Sensing" begin
        for FT in [Float32, Float64]
            can = ClimaCache.HyperspectralMLCanopy{FT}();
            for rs in [CanopyRadiativeTransfer.MODIS_EVI(can),
                       CanopyRadiativeTransfer.MODIS_EVI2(can),
                       CanopyRadiativeTransfer.MODIS_LSWI(can),
                       CanopyRadiativeTransfer.MODIS_NDVI(can),
                       CanopyRadiativeTransfer.MODIS_NIRv(can),
                       CanopyRadiativeTransfer.OCO2_SIF759(can),
                       CanopyRadiativeTransfer.OCO2_SIF770(can),
                       CanopyRadiativeTransfer.OCO3_SIF759(can),
                       CanopyRadiativeTransfer.OCO3_SIF770(can),
                       CanopyRadiativeTransfer.TROPOMI_SIF683(can),
                       CanopyRadiativeTransfer.TROPOMI_SIF740(can),
                       CanopyRadiativeTransfer.TROPOMI_SIF747(can),
                       CanopyRadiativeTransfer.TROPOMI_SIF751(can)]
                @test typeof(rs) == FT;
            end;
        end;
    end;
end;
