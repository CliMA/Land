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
            for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                result = envir_diff!(FT(0.1), mod, can, hs, envir, sm, 1);
                recursive_FT_test(result, FT);
                recursive_NaN_test(result);
            end
        end

        # test the stomata solutions
        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                leaf_photo_from_envir!(mod, can, hs, envir, sm, 1);
                recursive_NaN_test(can);
            end
        end

        for (mod,can) in zip([mod_3, mod_3], [can_3, can_4])
            for sm in [esm_1, esm_2, esm_3, esm_4, osm_1, osm_2, osm_3, osm_4, osm_5]
                leaf_photo_from_envir!(mod, can, hs, envir, sm);
                recursive_NaN_test(can);
            end
        end
    end
end
