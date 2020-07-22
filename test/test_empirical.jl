# test the structures
println("\nTesting and benchmarking empirical formulations...")
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
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            _beta = FT(1)
            @btime empirical_gsw_from_model($esm_1, $leaf_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_2, $leaf_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_3, $leaf_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_4, $leaf_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_1, $can_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_2, $can_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_3, $can_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_4, $can_3, $envir, $_beta)
            @btime empirical_gsw_from_model($esm_1, $can_3, $envir, $_beta, 1)
            @btime empirical_gsw_from_model($esm_2, $can_3, $envir, $_beta, 1)
            @btime empirical_gsw_from_model($esm_3, $can_3, $envir, $_beta, 1)
            @btime empirical_gsw_from_model($esm_4, $can_3, $envir, $_beta, 1)
        end
    end
end
