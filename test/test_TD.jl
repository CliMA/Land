# FT and NaN test for the structs
println("Testing and benchmarking temperature dependence functions...");
@testset "Photosynthesis --- temperature dependence" begin
    for FT in [Float32, Float64]
        c3_set = C3CLM(FT);
        c4_set = C4CLM(FT);
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        envir  = AirLayer{FT}();
        T      = rand(FT) + 298;

        leaf_temperature_dependence!(c3_set, leaf_3, envir, T);
        recursive_FT_test(leaf_3, FT);
        recursive_NaN_test(leaf_3);

        leaf_temperature_dependence!(c4_set, leaf_4, envir, T);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_4);

        if benchmarking
            @btime leaf_temperature_dependence!($c3_set, $leaf_3, $envir, $T);
            @btime leaf_temperature_dependence!($c4_set, $leaf_4, $envir, $T);
        end
    end
end
