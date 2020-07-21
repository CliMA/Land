# FT and NaN test for the structs
println("Testing and benchmarking ETR functions...")
@testset "Photosynthesis --- ETR functions" begin
    for FT in [Float32, Float64]
        c3_set = C3CLM(FT);
        c4_set = C4CLM(FT);
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        envir  = AirLayer{FT}();

        leaf_ETR!(c3_set, leaf_3);
        leaf_ETR!(c4_set, leaf_4);
        recursive_FT_test(leaf_3, FT);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_3);
        recursive_NaN_test(leaf_4);

        if benchmarking
            @btime leaf_ETR!($c3_set, $leaf_3);
            @btime leaf_ETR!($c4_set, $leaf_4);
        end
    end
end
