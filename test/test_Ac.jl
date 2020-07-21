# FT and NaN test for the structs
println("Testing and benchmarking rubisco_limited* functions...")
@testset "Photosynthesis --- rubisco_limited* functions" begin
    for FT in [Float32, Float64]
        c3_set = C3CLM(FT);
        c4_set = C4CLM(FT);
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        envir  = AirLayer{FT}();

        rubisco_limited_rate!(c3_set, leaf_3);
        recursive_FT_test(leaf_3, FT);
        recursive_NaN_test(leaf_3);

        rubisco_limited_rate!(c4_set, leaf_4);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_4);

        rubisco_limited_rate_glc!(c3_set, leaf_3, envir);
        recursive_FT_test(leaf_3, FT);
        recursive_NaN_test(leaf_3);

        if benchmarking
            @btime rubisco_limited_rate!($c3_set, $leaf_3);
            @btime rubisco_limited_rate!($c4_set, $leaf_4);
            @btime rubisco_limited_rate_glc!($c3_set, $leaf_3, $envir);
        end
    end
end
