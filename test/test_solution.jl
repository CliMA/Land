# test the structures
println("\nTesting and benchmarking stomatal models...")
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
                    recursive_FT_test(result, FT);
                    recursive_NaN_test(result);
                end
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                result = envir_diff!(FT(0.1), mod, can, hs, envir, sm, 1);
                recursive_FT_test(result, FT);
                recursive_NaN_test(result);
            end
        end

        # test the stomata solutions
        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4]
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}(), 1);
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}(), 1);
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}(), 1);
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}(), 1);
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}(), 1);
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}(), 1);
                recursive_NaN_test(can);
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                leaf_photo_from_envir!(mod, can, hs, envir, sm, 1);
                recursive_NaN_test(can);
            end
        end

        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4]
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}());
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}());
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}());
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}());
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}());
                recursive_NaN_test(can);
                leaf_photo_from_envir!(mod, can, hs, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}());
                recursive_NaN_test(can);
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                leaf_photo_from_envir!(mod, can, hs, envir, sm);
                recursive_NaN_test(can);
            end
        end

        # benchmarking the stomata solutions
        if benchmarking
            # test the stomata solutions
            println("Benchmarking the function for individual leaf...");
            for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
                for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                    @btime leaf_photo_from_envir!($mod, $can, $hs, $envir, $sm, 1);
                end
            end

            println("Benchmarking the function for whole canopy layer...");
            for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
                for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                    @btime leaf_photo_from_envir!($mod, $can, $hs, $envir, $sm);
                end
            end
        end
    end
end
