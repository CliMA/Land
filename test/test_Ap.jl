# FT and NaN test for the structs
println("Testing and benchmarking product_limited* functions...")
@testset "Photosynthesis --- product_limited* functions" begin
    for FT in [Float32, Float64]
        c3_set = C3CLM(FT);
        c4_set = C4CLM(FT);
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        envir  = AirLayer{FT}();

        product_limited_rate!(c3_set, leaf_3);
        recursive_FT_test(leaf_3, FT);
        recursive_NaN_test(leaf_3);

        product_limited_rate!(c4_set, leaf_4);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_4);

        product_limited_rate_glc!(c4_set, leaf_4, envir);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_4);

        if benchmarking
            @btime product_limited_rate!($c3_set, $leaf_3);
            @btime product_limited_rate!($c4_set, $leaf_4);
            @btime product_limited_rate_glc!($c4_set, $leaf_4, $envir);
        end
    end
end
