# test the structures
println("\nTesting and benchmarking gas exchange functions...")
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
        recursive_NaN_test(can_3);
        recursive_NaN_test(can_4);

        update_leaf_AK!(mod_3, can_3, hs, envir);
        update_leaf_AK!(mod_4, can_4, hs, envir);
        recursive_NaN_test(can_3);
        recursive_NaN_test(can_4);

        update_leaf_from_glc!(mod_3, can_3, envir, 1, FT(0.1));
        update_leaf_from_glc!(mod_4, can_4, envir, 1, FT(0.1));
        recursive_NaN_test(can_3);
        recursive_NaN_test(can_4);

        update_leaf_from_gsw!(mod_3, can_3, envir, 1, FT(0.05));
        update_leaf_from_gsw!(mod_4, can_4, envir, 1, FT(0.05));
        recursive_NaN_test(can_3);
        recursive_NaN_test(can_4);

        can_3.g_sw[2] = 0;
        can_4.g_sw[2] = 0;
        leaf_gsw_control!(mod_3, can_3, envir, 2);
        leaf_gsw_control!(mod_4, can_4, envir, 2);
        recursive_NaN_test(can_3);
        recursive_NaN_test(can_4);

        if benchmarking
            _glc = FT(0.1);
            _gsw = FT(0.1);
            @btime update_leaf_TP!($mod_3, $can_3, $hs, $envir);
            @btime update_leaf_TP!($mod_4, $can_4, $hs, $envir);
            @btime update_leaf_AK!($mod_3, $can_3, $hs, $envir);
            @btime update_leaf_AK!($mod_4, $can_4, $hs, $envir);
            @btime update_leaf_from_glc!($mod_3, $can_3, $envir, 1, $_glc);
            @btime update_leaf_from_glc!($mod_4, $can_4, $envir, 1, $_glc);
            @btime update_leaf_from_gsw!($mod_3, $can_3, $envir, 1, $_gsw);
            @btime update_leaf_from_gsw!($mod_4, $can_4, $envir, 1, $_gsw);
            can_3.g_sw[2] = 0;
            can_4.g_sw[2] = 0;
            @btime leaf_gsw_control!($mod_3, $can_3, $envir, 2);
            @btime leaf_gsw_control!($mod_4, $can_4, $envir, 2);
        end
    end
end
