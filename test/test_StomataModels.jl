# test the structures
println("\nTesting FT and NaN for the structures...")
@testset "StomtaModels --- structures" begin
    for FT in [Float32, Float64]
        can_3 = CanopyLayer{FT}(n_leaf=2);
        esm_1 = ESMBallBerry{FT}();
        esm_2 = ESMGentine{FT}();
        esm_3 = ESMLeuning{FT}();
        esm_4 = ESMMedlyn{FT}();
        osm_1 = OSMEller{FT}();
        osm_2 = OSMSperry{FT}();
        osm_3 = OSMWang{FT}();
        osm_4 = OSMWAP{FT}();
        osm_5 = OSMWAPMod{FT}();

        # test the structures
        for result in [ can_3, esm_1, esm_2, esm_3, esm_4, osm_4, osm_5 ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end



# test the structures
println("\nTesting gas exchange functions...")
@testset "StomtaModels --- gas exchange functions" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        mod_3  = C3CLM(FT);
        mod_4  = C4CLM(FT);
        rand_T = rand(FT) + 298;
        can_3  = CanopyLayer{FT}(n_leaf=2);
        can_4  = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the refresh functions
        update_leaf_TP!(mod_3, can_3, hs, envir);
        update_leaf_TP!(mod_4, can_4, hs, envir);
        @test NaN_test(can_3);
        @test NaN_test(can_4);

        update_leaf_AK!(mod_3, can_3, hs, envir);
        update_leaf_AK!(mod_4, can_4, hs, envir);
        @test NaN_test(can_3);
        @test NaN_test(can_4);

        update_leaf_from_glc!(mod_3, can_3, envir, 1, FT(0.1));
        update_leaf_from_glc!(mod_4, can_4, envir, 1, FT(0.1));
        @test NaN_test(can_3);
        @test NaN_test(can_4);

        update_leaf_from_gsw!(mod_3, can_3, envir, 1, FT(0.05));
        update_leaf_from_gsw!(mod_4, can_4, envir, 1, FT(0.05));
        @test NaN_test(can_3);
        @test NaN_test(can_4);

        can_3.g_sw[2] = 0;
        can_4.g_sw[2] = 0;
        leaf_gsw_control!(mod_3, can_3, envir, 2);
        leaf_gsw_control!(mod_4, can_4, envir, 2);
        @test NaN_test(can_3);
        @test NaN_test(can_4);
    end
end



# test the structures
println("\nTesting empirical formulations...")
@testset "StomtaModels --- empirical formulations" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        mod_3  = C3CLM(FT);
        mod_4  = C4CLM(FT);
        rand_T = rand(FT) + 298;
        can_3  = CanopyLayer{FT}(n_leaf=2);
        can_4  = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the empirical model formulations
        for result in [ empirical_gsw_from_model(esm_1, leaf_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_2, leaf_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_3, leaf_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_4, leaf_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_1, can_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_2, can_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_3, can_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_4, can_3, envir, FT(1)),
                        empirical_gsw_from_model(esm_1, can_3, envir, FT(1), 1),
                        empirical_gsw_from_model(esm_2, can_3, envir, FT(1), 1),
                        empirical_gsw_from_model(esm_3, can_3, envir, FT(1), 1),
                        empirical_gsw_from_model(esm_4, can_3, envir, FT(1), 1) ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end



# test the structures
println("\nTesting stomatal models...")
@testset "StomtaModels --- stomatal models" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        mod_3  = C3CLM(FT);
        mod_4  = C4CLM(FT);
        rand_T = rand(FT) + 298;
        can_3  = CanopyLayer{FT}(n_leaf=2);
        can_4  = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the solution functions
        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4]
                for result in [ envir_diff!(FT(0.1), mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}(), 1),
                                envir_diff!(FT(0.1), mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}(), 1),
                                envir_diff!(FT(0.1), mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}(), 1),
                                envir_diff!(FT(0.1), mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}(), 1),
                                envir_diff!(FT(0.1), mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}(), 1),
                                envir_diff!(FT(0.1), mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}(), 1),
                                ]
                    @test FT_test(result, FT);
                    @test NaN_test(result);
                end
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                result = envir_diff!(FT(0.1), mod, can, hs, envir, sm, 1);
                @test FT_test(result, FT);
                @test NaN_test(result);
            end
        end

        # test the stomata solutions
        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4]
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}(), 1);
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}(), 1);
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}(), 1);
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}(), 1);
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}(), 1);
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}(), 1);
                @test NaN_test(can);
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                leaf_photo_from_envir!(mod, can, hs, envir, sm, 1);
                @test NaN_test(can);
            end
        end

        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4]
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}());
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}());
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}());
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}());
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}());
                @test NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}());
                @test NaN_test(can);
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                leaf_photo_from_envir!(mod, can, hs, envir, sm);
                @test NaN_test(can);
            end
        end
    end
end
