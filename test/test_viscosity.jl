# test and benchmark laten heat functions
println("\nTesting and Benchmarking viscosity* functions...");
@testset "Testing + Benchmarking --- viscosity*" begin
    for FT in [Float32, Float64]
        rand_T  = rand(FT) + 298
        rand_Tl = rand(FT,10) .+ 298
        for result in [ viscosity(rand_T),
                        viscosity.(rand_Tl),
                        relative_viscosity(rand_T),
                        relative_viscosity.(rand_Tl) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        @test relative_viscosity(FT(298.15)) ≈ 1;

        if benchmarking
            @show FT;
            @btime viscosity($rand_T);
            @btime viscosity.($rand_Tl);
            @btime relative_viscosity($rand_T);
            @btime relative_viscosity.($rand_Tl);
        end
    end
end
