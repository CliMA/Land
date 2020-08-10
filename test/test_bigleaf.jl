# Test big leaf canopy model
println("\nTesting and benchmarking the big leaf model...");
@testset "CanopyRadiation --- big leaf model" begin
    for FT in [Float32, Float64]
        for result in [ big_leaf_partition(FT(3.0), FT(30.0), FT(1000.0)),
                        big_leaf_partition(FT(2.0), FT(30.0), FT(1000.0)),
                        big_leaf_partition(FT(1.0), FT(30.0), FT(1000.0)) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            _lai  = FT(2.0);
            _zen  = FT(30.0);
            _rall = FT(1000.0);
            @btime big_leaf_partition($_lai, $_zen, $_rall);
        end
    end
end
