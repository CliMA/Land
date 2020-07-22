# FT and NaN test
println("\nTesting the structures...")
@testset "Hydraulics --- structures" begin
    for FT in [Float32, Float64]
        leaf  = LeafHydraulics{FT}();
        root  = RootHydraulics{FT}();
        stem  = StemHydraulics{FT}();
        grass = create_grass_like_hs(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
        palm  = create_palm_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        tree  = create_tree_like_hs(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
        treet = TreeSimple{FT}();
        _vc1  = WeibullSingle{FT}();
        _vc2  = WeibullDual{FT}();
        _sh1  = BrooksCorey{FT}();
        _sh2  = VanGenuchten{FT}();

        # Test the struct
        for data_set in [ leaf, root, stem, grass, palm, tree, treet, _vc1, _vc2, _sh1, _sh2 ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set)
        end
    end
end
