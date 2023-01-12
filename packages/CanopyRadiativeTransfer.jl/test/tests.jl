@testset verbose = true "CanopyRadiativeTransfer Test" begin
    @testset "Leaf Inclination Angles" begin
        for FT in [Float32, Float64]
            can1 = EmeraldNamespace.HyperspectralMLCanopy{FT}();
            can2 = EmeraldNamespace.BroadbandSLCanopy{FT}();
            CanopyRadiativeTransfer.inclination_angles!(can1, can1.LIDF);
            CanopyRadiativeTransfer.inclination_angles!(can2, can2.LIDF);
            @test true;
        end;
    end;

    @testset "Clumping Index" begin
        for FT in [Float32, Float64]
            can = EmeraldNamespace.HyperspectralMLCanopy{FT}();
            angles = EmeraldNamespace.SunSensorGeometry{FT}();
            CanopyRadiativeTransfer.clumping_index!(can, angles);
            @test true;
        end;
    end;

    @testset "Canopy RT" begin
        for FT in [Float32, Float64]
            hcan = EmeraldNamespace.HyperspectralMLCanopy{FT}();
            bcan = EmeraldNamespace.BroadbandSLCanopy{FT}();
            angles = EmeraldNamespace.SunSensorGeometry{FT}();
            hleaf = EmeraldNamespace.Leaves2D{FT}();
            bleaf = EmeraldNamespace.Leaves1D{FT}();
            leaves = [deepcopy(hleaf) for i in 1:20];
            hsoil = EmeraldNamespace.Soil{FT}();
            bsoil = EmeraldNamespace.Soil{FT}(ZS = FT[0,-1], ALBEDO = EmeraldNamespace.BroadbandSoilAlbedo{FT}());
            hrad = EmeraldNamespace.HyperspectralRadiation{FT}();
            brad = EmeraldNamespace.BroadbandRadiation{FT}();
            spac = EmeraldNamespace.MonoMLTreeSPAC{FT}();
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
            can = EmeraldNamespace.HyperspectralMLCanopy{FT}();
            for var in [CanopyRadiativeTransfer.MODIS_EVI(can),
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
                @test typeof(var) == FT;
            end;
        end;
    end;
end;
