using ClimaCache
using StomataModels
using Test


@testset verbose = true "StomataModels Test" begin
    @testset "Beta Function" begin
        for FT in [Float32, Float64]
            f_1(x) = x;
            f_2(x) = min(1, max(0, (x-0.1) / 0.3));
            xvc = ClimaCache.WeibullVC{FT}(2,5);
            svc = ClimaCache.VanGenuchten{FT}("Loam");
            swc = FT(0.3)
            @test 0 < StomataModels.β_factor(f_1, xvc, FT(-1)) < 1;
            @test 0 < StomataModels.β_factor(f_1, svc, FT(-1)) < 1;
            @test 0 < StomataModels.β_factor(f_2, swc) < 1;
        end;
    end;
end;
