using PlantHydraulics
using Test


@testset verbose = true "PlantHydraulics Test" begin
    @testset "Vulnerability" begin
        for FT in [Float32, Float64]
            vc = LogisticVC{FT}(2, 2);
            k1 = xylem_k_ratio(vc, FT(-2.0));
            @test true;
            k2 = xylem_k_ratio(vc, FT(-2.0), FT(0.95));
            @test true;
        end;
    end;
end;
