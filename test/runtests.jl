using CanopyRadiativeTransfer
using ClimaCache
using Test


@testset verbose = true "CanopyRadiativeTransfer Test" begin
    @testset "Leaf Inclination Angles" begin
        for FT in [Float32, Float64]
            can1 = ClimaCache.HyperspectralMLCanopy{FT}();
            can2 = ClimaCache.BroadbandSLCanopy{FT}();
            CanopyRadiativeTransfer.inclination_angles!(can1, can1.LIDF, FT(0), FT(0));
            CanopyRadiativeTransfer.inclination_angles!(can2, can2.LIDF, FT(0), FT(0));
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

    @testset "Canopy RT" begin
        for FT in [Float32, Float64]
            hcan = ClimaCache.HyperspectralMLCanopy{FT}();
            bcan = ClimaCache.BroadbandSLCanopy{FT}();
            angles = ClimaCache.SunSensorGeometry{FT}();
            hleaf = ClimaCache.Leaves2D{FT}("C3");
            bleaf = ClimaCache.Leaves1D{FT}("C3");
            leaves = [deepcopy(hleaf) for i in 1:20];
            hsoil = ClimaCache.Soil{FT}(FT[0,-1]);
            bsoil = ClimaCache.Soil{FT}(FT[0,-1], true);
            hrad = ClimaCache.HyperspectralRadiation{FT}();
            brad = ClimaCache.BroadbandRadiation{FT}();
            spac = ClimaCache.MonoMLTreeSPAC{FT}("C3");
            CanopyRadiativeTransfer.canopy_optical_properties!(hcan, angles);
            @test true;
            CanopyRadiativeTransfer.canopy_optical_properties!(hcan, leaves, hsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_optical_properties!(hcan, leaves, bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(hcan, leaves, hrad, hsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(hcan, leaves, hrad, bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(hcan, leaves, FT(100), hsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(hcan, leaves, FT(100), bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(bcan, bleaf, brad, bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(bcan, bleaf, FT(100), bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_fluorescence!(hcan, leaves);
            @test true;
            CanopyRadiativeTransfer.soil_albedo!(hcan, hsoil);
            @test true;
            CanopyRadiativeTransfer.soil_albedo!(hcan, bsoil);
            @test true;
            CanopyRadiativeTransfer.canopy_radiation!(spac);
            @test true;
            CanopyRadiativeTransfer.canopy_fluorescence!(spac);
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
