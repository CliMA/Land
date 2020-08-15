# test and benchmark laten heat functions
println("\nTesting and Benchmarking capillary_pressure functions...");
@testset "Testing + Benchmarking --- capillary_pressure" begin
    for FT in [Float32, Float64]
        rand_r = (rand(FT) + 20) * FT(1e-6);
        rand_T = rand(FT) + 298;
        rand_α = rand(FT) * 50;
        for result in [ capillary_pressure(rand_r, rand_T),
                        capillary_pressure(rand_r, rand_T, rand_α) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @show FT;
            @btime capillary_pressure($rand_r, $rand_T);
            @btime capillary_pressure($rand_r, $rand_T, $rand_α);
        end
    end
end
