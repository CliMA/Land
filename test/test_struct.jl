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
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end
end
