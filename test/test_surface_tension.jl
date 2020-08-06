# test and benchmark laten heat functions
println("\nTesting and Benchmarking surface_tension* functions...");
@testset "Testing + Benchmarking --- surface_tension*" begin
    for FT in [Float32, Float64]
        rand_T  = rand(FT) + 298
        rand_Tl = rand(FT,10) .+ 298
        for result in [ surface_tension(rand_T),
                        surface_tension.(rand_Tl),
                        relative_surface_tension(rand_T),
                        relative_surface_tension.(rand_Tl) ]
            recursive_FT_test(result, FT);
            recursive_NaN_test(result);
        end

        @test relative_surface_tension(FT(298.15)) ≈ 1;

        if benchmarking
            @show FT;
            @btime surface_tension($rand_T);
            @btime surface_tension.($rand_Tl);
            @btime relative_surface_tension($rand_T);
            @btime relative_surface_tension.($rand_Tl);
        end
    end
end
