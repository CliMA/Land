# test and benchmark laten heat functions
println("\nTesting and Benchmarking relative_diffusive_coefficient functions...");
@testset "Testing + Benchmarking --- relative_diffusive_coefficient" begin
    for FT in [Float32, Float64]
        rand_T = rand(FT) + 298
        for result in [ relative_diffusive_coefficient(rand_T) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        @test relative_diffusive_coefficient(FT(298.15)) ≈ 1;

        if benchmarking
            @show FT;
            @btime relative_diffusive_coefficient($rand_T);
        end
    end
end
