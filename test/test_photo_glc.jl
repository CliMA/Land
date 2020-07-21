# FT and NaN test for the structs
println("Testing and benchmarking leaf_photo_from_glc! functions...")
@testset "Photosynthesis --- leaf_photo_from_glc! functions" begin
    for FT in [Float32, Float64]
        c3_set = C3CLM(FT);
        c4_set = C4CLM(FT);
        leaf_3 = Leaf{FT}();
        leaf_4 = Leaf{FT}();
        envir  = AirLayer{FT}();
        T      = rand(FT) + 298;
        glc    = FT(0.1);

        leaf_temperature_dependence!(c3_set, leaf_3, envir, T);
        leaf_temperature_dependence!(c4_set, leaf_4, envir, T);

        leaf_photo_from_glc!(c3_set, leaf_3, envir, glc);
        recursive_FT_test(leaf_3, FT);
        recursive_NaN_test(leaf_3);

        leaf_photo_from_glc!(c4_set, leaf_4, envir, glc);
        recursive_FT_test(leaf_4, FT);
        recursive_NaN_test(leaf_4);

        if benchmarking
            @btime leaf_photo_from_glc!($c3_set, $leaf_3, $envir, $glc);
            @btime leaf_photo_from_glc!($c4_set, $leaf_4, $envir, $glc);
        end
    end
end
