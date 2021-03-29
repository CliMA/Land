# FT and NaN test for the structs
@testset "Photosynthesis --- structures" begin
    for FT in [Float32, Float64]
        for data_set in [ KcTDBernacchi(FT),
                          VpmaxTDBoyd(FT),
                          C3CLM(FT),
                          C4CLM(FT),
                          AirLayer{FT}(),
                          Leaf{FT}(),
                          FluorescenceVanDerTol(FT),
                          FluorescenceVanDerTolDrought(FT),
                          KoTDBernacchi(FT),
                          RespirationTDBernacchi(FT),
                          VcmaxTDBernacchi(FT),
                          VomaxTDBernacchi(FT),
                          ΓStarTDBernacchi(FT),
                          KpepTDBoyd(FT),
                          JmaxTDLeuning(FT),
                          VcmaxTDLeuning(FT),
                          JmaxTDBernacchi(FT),
                          VtoRCollatz(FT),
                          C3Bernacchi(FT),
                          Q10TDAngiosperm(FT),
                          Q10TDGymnosperm(FT) ]
            @test FT_test(data_set, FT);
            @test NaN_test(data_set);
        end
    end
end




# FT and NaN test for the structs
println();
@testset "Photosynthesis --- functions" begin
    for FT in [Float32, Float64]
        c3_set   = C3CLM(FT);
        c4_set   = C4CLM(FT);
        leaf_3   = Leaf{FT}();
        leaf_4   = Leaf{FT}();
        td_q10   = Q10TD{FT}(1, 273.15, 1.7);
        envir    = AirLayer{FT}();
        fluo_set = c3_set.Flu;
        T        = rand(FT) + 298;
        glc      = FT(0.1);
        p_i      = rand(FT) + 20;

        # temperature corrections
        photo_TD_from_set(td_q10, T);
        leaf_temperature_dependence!(c3_set, leaf_3, envir, T);
        leaf_temperature_dependence!(c4_set, leaf_4, envir, T);

        # rubisco limited rates
        rubisco_limited_rate!(c3_set, leaf_3);
        rubisco_limited_rate!(c4_set, leaf_4);
        @test NaN_test(leaf_3);
        @test NaN_test(leaf_4);
        rubisco_limited_rate!(c3_set, leaf_3, envir);
        @test NaN_test(leaf_3);

        # light limited rates
        leaf_ETR!(c3_set, leaf_3);
        leaf_ETR!(c4_set, leaf_4);
        light_limited_rate!(c3_set, leaf_3);
        light_limited_rate!(c4_set, leaf_4);
        @test NaN_test(leaf_3);
        @test NaN_test(leaf_4);
        light_limited_rate!(c3_set, leaf_3, envir);
        @test NaN_test(leaf_3);

        # product limited rates
        product_limited_rate!(c3_set, leaf_3);
        product_limited_rate!(c4_set, leaf_4);
        @test NaN_test(leaf_3);
        @test NaN_test(leaf_4);
        product_limited_rate!(c4_set, leaf_4, envir);
        @test NaN_test(leaf_4);

        # fluorescence
        leaf_photosynthesis!(c3_set, leaf_3, envir, PCO₂Mode(), FT(2));
        leaf_fluorescence!(fluo_set, leaf_3);
        leaf_photosynthesis!(c3_set, leaf_3, envir, GCO₂Mode());
        leaf_fluorescence!(fluo_set, leaf_3);
        @test NaN_test(leaf_3);

        # leaf photo from glc
        leaf_photosynthesis!(c3_set, leaf_3, envir, GCO₂Mode(), glc);
        leaf_photosynthesis!(c4_set, leaf_4, envir, GCO₂Mode(), glc);
        @test NaN_test(leaf_3);
        @test NaN_test(leaf_4);

        # leaf photo from p_i
        leaf_photosynthesis!(c3_set, leaf_3, envir, PCO₂Mode(), p_i);
        leaf_photosynthesis!(c4_set, leaf_4, envir, PCO₂Mode(), p_i);
        @test NaN_test(leaf_3);
        @test NaN_test(leaf_4);
    end
end
