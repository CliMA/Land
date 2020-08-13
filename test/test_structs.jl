# FT and NaN tests
println("\nTesting the structures...");
@testset "CanopyRadiation --- FT and NaN test" begin
    for FT in [Float32, Float64]
        for data_set in initialize_rt_module(FT)
            recursive_FT_test(data_set, FT);
            recursive_NaN_test(data_set);
        end
    end
end
