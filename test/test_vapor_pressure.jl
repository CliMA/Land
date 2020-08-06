# test and benchmark laten heat functions
println("\nTesting and Benchmarking saturation_vapor_pressure* functions...");
@testset "Testing + Benchmarking --- saturation_vapor_pressure*" begin
    for FT in [Float32, Float64]
        rand_T  = rand(FT) + 298
        rand_Tl = rand(FT,10) .+ 298
        for result in [ saturation_vapor_pressure(rand_T),
                        saturation_vapor_pressure.(rand_Tl),
                        saturation_vapor_pressure_slope(rand_T),
                        saturation_vapor_pressure_slope.(rand_Tl) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        if benchmarking
            @show FT;
            @btime saturation_vapor_pressure($rand_T);
            @btime saturation_vapor_pressure.($rand_Tl);
            @btime saturation_vapor_pressure_slope($rand_T);
            @btime saturation_vapor_pressure_slope.($rand_Tl);
        end
    end
end
