# FT and NaN test
println("\nTesting and benchmarking the plant-level functions...")
@testset "Hydraulics --- plant-level" begin
    for FT in [Float32, Float64]
        grass = create_grass_like_hs(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
        palm  = create_palm_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        tree  = create_tree_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        treet = TreeSimple{FT}();

        # test the tree_e_crit function
        for result in [ tree_e_crit(grass),
                        tree_e_crit(palm),
                        tree_e_crit(tree),
                        tree_e_crit(treet) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime tree_e_crit($grass);
            @btime tree_e_crit($palm);
            @btime tree_e_crit($tree);
            @btime tree_e_crit($treet);
        end
    end
end
