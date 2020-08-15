# test and benchmark laten heat functions
println("\nTesting and Benchmarking latent_heat_vapor functions...");
@testset "Testing + Benchmarking --- latent_heat_vapor" begin
    for FT in [Float32, Float64]
        rand_T = rand(FT) + 298;
        for result in [ latent_heat_vapor(rand_T) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @show FT;
            @btime latent_heat_vapor($rand_T);
        end
    end
end
