# Test big leaf canopy model
println("\nTesting the big leaf model...");
@testset "CanopyRadiation --- big leaf model" begin
    for FT in [Float32, Float64]
        for result in [ big_leaf_partition(FT(3.0), FT(30.0), FT(1000.0)),
                        big_leaf_partition(FT(2.0), FT(30.0), FT(1000.0)),
                        big_leaf_partition(FT(1.0), FT(30.0), FT(1000.0)) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end
end
