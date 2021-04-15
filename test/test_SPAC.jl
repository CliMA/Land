# test the structs
@info "Testing the FT and NaN of the structs...";
@testset "FT and NaN --- Types" begin
    for FT in [Float32, Float64]
        node = SPACSimple{FT}();
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test and benchmark the big_leaf_partition!
println();
@info "Testing the big_leaf_partition!...";
@testset "SoilPlantAirContinuum --- big_leaf_partition!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        zenith = FT(30);
        r_all  = FT(1000);

        big_leaf_partition!(node, zenith, r_all);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test the gain_risk_map
println();
@info "Testing the gain_risk_map Function...";
@testset "SoilPlantAirContinuum --- gain_risk_map" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);

        big_leaf_partition!(node, zenith, r_all);
        mat = gain_risk_map(node, photo);
        @test FT_test(mat, FT);
        @test NaN_test(mat);
    end
end




# test and benchmark the leaf_gas_exchange_nonopt!
println();
@info "Testing the leaf_gas_exchange_nonopt! Functions...";
@testset "SoilPlantAirContinuum --- leaf_gas_exchange_nonopt!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);
        flow   = FT(4);
        f_sl   = FT(2.5);
        f_sh   = FT(1.5);

        big_leaf_partition!(node, zenith, r_all);
        leaf_gas_exchange_nonopt!(node, photo, flow);
        @test FT_test(node, FT);
        @test NaN_test(node);
        leaf_gas_exchange_nonopt!(node, photo, f_sl, f_sh);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test and benchmark the leaf_gas_exchange!
println();
@info "Testing the leaf_gas_exchange! Functions...";
@testset "SoilPlantAirContinuum --- leaf_gas_exchange!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);
        flow   = FT(4);
        f_sl   = FT(2.5);
        f_sh   = FT(1.5);

        big_leaf_partition!(node, zenith, r_all);
        leaf_gas_exchange!(node, photo, flow);
        @test FT_test(node, FT);
        @test NaN_test(node);
        leaf_gas_exchange!(node, photo, f_sl, f_sh);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test and benchmark the leaf_temperature*
println();
@info "Testing the leaf_temperature* Functions...";
@testset "SoilPlantAirContinuum --- leaf_temperature*" begin
    for FT in [Float32, Float64]
        node = SPACSimple{FT}();
        rad  = FT(300);
        flow = FT(4);

        for result in [ leaf_temperature(node, rad, flow),
                        leaf_temperature_shaded(node, rad, flow),
                        leaf_temperature_sunlit(node, rad, flow) ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end




# test and benchmark the optimize_flows!
println();
@info "Testing the optimize_flows! Functions...";
@testset "SoilPlantAirContinuum --- optimize_flows!" begin
    for FT in [Float32, Float64]
        node   = SPACSimple{FT}();
        photo  = C3CLM(FT);
        zenith = FT(30);
        r_all  = FT(1000);

        big_leaf_partition!(node, zenith, r_all);
        optimize_flows!(node, photo);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test and benchmark the atmosheric* functions
println();
@info "Testing the atmosheric* Functions...";
@testset "SoilPlantAirContinuum --- atmosheric*" begin
    for FT in [Float32, Float64]
        h = FT(1000);

        for result in [ atmospheric_pressure(h),
                        atmospheric_pressure_ratio(h),
                        ppm_to_Pa(h) ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end




# test and benchmark the zenith_angle
println();
@info "Testing the zenith_angle Functions...";
@testset "SoilPlantAirContinuum --- zenith_angle" begin
    for FT in [Float32, Float64]
        latd = FT(10);
        decd = FT(10);
        lhad = FT(10);
        day  = FT(100)
        hour = FT(13)
        minu = FT(30)

        for result in [ zenith_angle(latd, decd, lhad),
                        zenith_angle(latd, day, hour),
                        zenith_angle(latd, day, hour, minu) ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end




# test and benchmark the annual_profit
println();
@info "Testing the annual_profit Functions...";
@testset "SoilPlantAirContinuum --- annual_profit" begin
    arti = artifact"2020_leaf_invest_weather" * "/gs_sample.csv";
    weat = read_csv(arti);
    for FT in [Float32, Float64]
        node    = SPACSimple{FT}();
        photo   = C3CLM(FT);
        weatmat = Matrix{FT}(weat);

        gscp = annual_profit(node, photo, weatmat);
        @test FT_test(gscp, FT);
        @test NaN_test(gscp);
    end
end




# test and benchmark the annual_simulation!
println();
@info "Testing annual_simulation! Functions...";
@testset "SoilPlantAirContinuum --- annual_simulation!" begin
    arti = artifact"2020_leaf_invest_weather" * "/gs_sample.csv";
    weat = read_csv(arti);
    for FT in [Float32, Float64]
        node  = SPACSimple{FT}();
        photo = C3CLM(FT);
        df    = create_dataframe(FT, weat);

        annual_simulation!(node, photo, weat, df);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test and benchmark the leaf_allocation!
println();
@info "Testing the leaf_allocation! Functions...";
@testset "SoilPlantAirContinuum --- leaf_allocation!" begin
    for FT in [Float32, Float64]
        node  = SPACSimple{FT}();
        photo = C3CLM(FT);
        laba  = FT(1000);
        vmax  = FT(80);

        leaf_allocation!(node, laba);
        @test FT_test(node, FT);
        @test NaN_test(node);
        leaf_allocation!(node, photo, vmax);
        @test FT_test(node, FT);
        @test NaN_test(node);
        leaf_allocation!(node, photo, laba, vmax);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test and benchmark the optimize_leaf
println();
@info "Testing the optimize_leaf! Functions...";
@testset "SoilPlantAirContinuum --- optimize_leaf!" begin
    arti = artifact"2020_leaf_invest_weather" * "/gs_sample.csv";
    weat = read_csv(arti);
    for FT in [Float32, Float64]
        node    = SPACSimple{FT}();
        photo   = C3CLM(FT);
        weatmat = Matrix{FT}(weat);

        optimize_leaf!(node, photo, weatmat);
        @test FT_test(node, FT);
        @test NaN_test(node);
    end
end




# test the function to vary SPACSimple
println();
@info "Testing the vary_spac! Functions...";
@testset "SoilPlantAirContinuum --- vary_spac!" begin
    arti = artifact"2020_leaf_invest_weather" * "/gs_sample.csv";
    weat = read_csv(arti);
    facs = ["kl", "kw", "wb", "wc", "wk",
            "cc", "cv", "gm",
            "ga", "sd",
            "ta", "rh", "ca"];
    for FT in [Float32, Float64]
        node = SPACSimple{FT}();
        for _fac in facs
            vary_spac!(node, weat, _fac, FT(1.5));
            @test true;
        end
    end
end




# test the function for mSCOPE version radiation
println();
@info "Testing the mSCOPE verion setup...";
@testset "SoilPlantAirContinuum --- vary_spac!" begin
    for FT in [Float32, Float64]
        node = SPACMono{FT}();
        initialize_spac_canopy!(node);
        layer_fluxes!(node);
        @test NaN_test(node);

        update_Cab!(node, FT(30));
        update_Kmax!(node, FT(1));
        update_LAI!(node, FT(3));
        update_VJR!(node, FT(0.5));
        update_VJRWW!(node, FT(50));
        update_Weibull!(node, FT(3));
        update_Weibull!(node, FT(3), FT(0.9));
        @test true;
    end
end
