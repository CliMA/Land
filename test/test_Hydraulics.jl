# FT and NaN test
@testset "Hydraulics --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        leaf = PH.LeafHydraulics{FT}();
        stem = PH.StemHydraulics{FT}();

        for data_set in [ leaf,
                          stem ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set)
        end

        _f_1 = FT(0.001)
        _f_2 = FT(1)
        for result in [ PH.xylem_p_from_flow(leaf, _f_1),
                        PH.xylem_p_from_flow(leaf, _f_2),
                        PH.leaf_xylem_risk(leaf, _f_1),
                        PH.leaf_xylem_risk(leaf, _f_2),
                        PH.leaf_e_crit(leaf) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end

end
