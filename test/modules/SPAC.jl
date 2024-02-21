@testset verbose = true "SoilPlantAirContinuum" begin
    @testset "Atmosheric pressure" begin
        for FT in [Float32, Float64]
            h = FT(1000);

            for result in [ SoilPlantAirContinuum.atmospheric_pressure(h),
                            SoilPlantAirContinuum.atmospheric_pressure_ratio(h),
                            SoilPlantAirContinuum.ppm_to_Pa(h) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Solar zenith_angle" begin
        for FT in [Float32, Float64]
            latd = FT(10);
            decd = FT(10);
            lhad = FT(10);
            day  = FT(100);
            hour = FT(13);
            minu = FT(30);

            for result in [ SoilPlantAirContinuum.zenith_angle(latd, decd, lhad),
                            SoilPlantAirContinuum.zenith_angle(latd, day, hour),
                            SoilPlantAirContinuum.zenith_angle(latd, day, hour, minu) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Vary SPAC parameters!" begin
        for FT in [Float32, Float64]
            node = SoilPlantAirContinuum.SPACMono{FT}();
            SoilPlantAirContinuum.initialize_spac_canopy!(node);
            SoilPlantAirContinuum.layer_fluxes!(node);
            SoilPlantAirContinuum.layer_fluxes!(node, FT(30));
            @test PkgUtility.NaN_test(node);

            SoilPlantAirContinuum.update_Cab!(node, FT(30));
            SoilPlantAirContinuum.update_Kmax!(node, FT(1));
            SoilPlantAirContinuum.update_LAI!(node, FT(3));
            SoilPlantAirContinuum.update_VJR!(node, FT(0.5));
            SoilPlantAirContinuum.update_VJRWW!(node, FT(50));
            SoilPlantAirContinuum.update_Weibull!(node, FT(3));
            SoilPlantAirContinuum.update_Weibull!(node, FT(3), FT(0.9));
            @test true;
        end;
    end;

    @testset "Switch among land spectra sets" begin
        for FT in [Float32, Float64]
            node = SoilPlantAirContinuum.SPACMono{FT}();
            SoilPlantAirContinuum.initialize_spac_canopy!(node);
            SoilPlantAirContinuum.layer_fluxes!(node);
            SoilPlantAirContinuum.layer_fluxes!(node, FT(30));
            ags = SoilPlantAirContinuum.A_GROSS(node);
            @test PkgUtility.NaN_test(ags);
            @test PkgUtility.NaN_test(node);
            @test true;
            node = SoilPlantAirContinuum.SPACMono{FT}(opti_file=CanopyLayers.LAND_2017);
            SoilPlantAirContinuum.initialize_spac_canopy!(node);
            SoilPlantAirContinuum.layer_fluxes!(node);
            SoilPlantAirContinuum.layer_fluxes!(node, FT(30));
            ags = SoilPlantAirContinuum.A_GROSS(node);
            @test PkgUtility.NaN_test(ags);
            @test PkgUtility.NaN_test(node);
            @test true;
        end;
    end;
end;
