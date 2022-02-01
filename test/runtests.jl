using PlantHydraulics
using Test


@testset verbose = true "PlantHydraulics Test" begin
    @testset "Vulnerability" begin
        for FT in [Float32, Float64]
            vc1 = LogisticVC{FT}(2, 2);
            kr1 = relative_hydraulic_conductance(vc1, FT(-2.0));
            @test true;

            vc2 = PowerVC{FT}(2, 2);
            kr2 = relative_hydraulic_conductance(vc2, FT(-2.0));
            @test true;

            vc3 = WeibullVC{FT}(2, 2);
            kr3 = relative_hydraulic_conductance(vc3, FT(-2.0));
            @test true;

            vs1 = ComplexVC{FT}([0.3,0.3,0.4], [vc1, vc1, vc1]);
            ks1 = relative_hydraulic_conductance(vs1, FT(-2.0));
            @test true;
            vs2 = ComplexVC{FT}([0.3,0.3,0.4], [vc2, vc2, vc2]);
            ks2 = relative_hydraulic_conductance(vs2, FT(-2.0));
            @test true;
            vs3 = ComplexVC{FT}([0.3,0.3,0.4], [vc3, vc3, vc3]);
            ks3 = relative_hydraulic_conductance(vs3, FT(-2.0));
            @test true;
            vs4 = ComplexVC{FT}([0.3,0.3,0.4], [vc1, vc2, vc3]);
            ks4 = relative_hydraulic_conductance(vs4, FT(-2.0));
            @test true;
        end;
    end;
end;
