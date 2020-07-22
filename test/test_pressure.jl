# FT and NaN test
println("\nTesting and benchmarking the pressure functions...")
@testset "Hydraulics --- pressure functions" begin
    for FT in [Float32, Float64]
        leaf  = LeafHydraulics{FT}();
        root  = RootHydraulics{FT}();
        stem  = StemHydraulics{FT}();
        grass = create_grass_like_hs(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
        palm  = create_palm_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        tree  = create_tree_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        treet = TreeSimple{FT}();

        # test the xylem_p_from_flow function
        _f_1 = FT(0.001)
        _f_2 = FT(1)
        for result in [ xylem_p_from_flow(leaf, _f_1),
                        xylem_p_from_flow(leaf, _f_2),
                        xylem_p_from_flow(root, _f_1),
                        xylem_p_from_flow(root, _f_2),
                        xylem_p_from_flow(stem, _f_1),
                        xylem_p_from_flow(stem, _f_2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime xylem_p_from_flow($leaf, $_f_1);
            @btime xylem_p_from_flow($root, $_f_1);
            @btime xylem_p_from_flow($stem, $_f_1);
        end

        # test the xylem_p_from_flow function for plants
        for result in [ xylem_p_from_flow(grass, _f_1),
                        xylem_p_from_flow(grass, _f_2),
                        xylem_p_from_flow(palm, _f_1),
                        xylem_p_from_flow(palm, _f_2),
                        xylem_p_from_flow(tree, _f_1),
                        xylem_p_from_flow(tree, _f_2),
                        xylem_p_from_flow(treet, _f_1),
                        xylem_p_from_flow(treet, _f_2) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime xylem_p_from_flow($grass, $_f_1);
            @btime xylem_p_from_flow($palm, $_f_1);
            @btime xylem_p_from_flow($tree, $_f_1);
            @btime xylem_p_from_flow($treet, $_f_1);
        end
    end
end
