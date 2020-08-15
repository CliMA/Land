# test and benchmark laten heat functions
println("\nTesting and Benchmarking saturation_vapor_pressure* functions...");
@testset "Testing + Benchmarking --- saturation_vapor_pressure*" begin
    for FT in [Float32, Float64]
        rand_T = rand(FT) + 298;
        rand_Ψ = rand(FT) - 3;
        for result in [ pressure_correction(rand_T, rand_Ψ),
                        saturation_vapor_pressure(rand_T),
                        saturation_vapor_pressure(rand_T, rand_Ψ),
                        saturation_vapor_pressure_slope(rand_T),
                        saturation_vapor_pressure_slope(rand_T, rand_Ψ) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end
    end
end
