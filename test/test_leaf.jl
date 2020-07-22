# FT and NaN test
println("\nTesting and benchmarking the leaf functions...")
@testset "Hydraulics --- leaf" begin
    for FT in [Float32, Float64]
        leaf  = LeafHydraulics{FT}();
        _f    = FT(0.01);

        # Test the struct
        _lr = leaf_xylem_risk(leaf, _f);
        _ec = leaf_e_crit(leaf);

        for result in [_lr, _ec]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @btime leaf_xylem_risk($leaf, $_f);
            @btime leaf_e_crit($leaf, $_ec);
        end
    end
end
