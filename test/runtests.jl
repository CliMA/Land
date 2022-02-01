using PlantHydraulics
using Test


@testset verbose = true "PlantHydraulics Test" begin
    @testset "Vulnerability" begin
        for FT in [Float32, Float64]
            vc = LogisticVC{FT}(2, 2);
            kr = relative_hydraulic_conductance(vc, FT(-2.0));
            @test true;

            vc = PowerVC{FT}(2, 2);
            kr = relative_hydraulic_conductance(vc, FT(-2.0));
            @test true;

            vc = WeibullVC{FT}(2, 2);
            kr = relative_hydraulic_conductance(vc, FT(-2.0));
            @test true;
        end;
    end;
end;
