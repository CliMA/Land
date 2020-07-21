# FT and NaN test for the structs
println("Testing and benchmarking Fluorescence function...")
@testset "Photosynthesis --- Fluorescence" begin
    for FT in [Float32, Float64]
        photo_set = C3CLM(FT);
        fluo_set  = photo_set.Flu;
        leaf      = Leaf{FT}();
        envir     = AirLayer{FT}();

        leaf_photo_from_glc!(photo_set, leaf, envir);
        leaf_fluorescence!(fluo_set, leaf);
        recursive_FT_test(leaf, FT);
        recursive_NaN_test(leaf);

        if benchmarking
            @btime leaf_fluorescence!($fluo_set, $leaf);
        end
    end
end
