using SoilHydraulics
using Test


@testset verbose = true "SoilHydraulics Test" begin
    @testset "Constructor" begin
        for FT in [Float32, Float64]
            bc = BrooksCorey{FT}(VanGenuchten{FT}("Loam"));
            @test true;
        end;
    end;

    @testset "Potential" begin
        for FT in [Float32, Float64]
            vg = VanGenuchten{FT}("Loam");
            bc = BrooksCorey{FT}(vg);
            @test soil_ψ_25(vg, FT(0.2)) < 0;
            @test soil_ψ_25(bc, FT(0.2)) < 0;
            @test soil_ψ_25(vg, FT(1.0)) == 0;
            @test soil_ψ_25(bc, FT(1.0)) == 0;
        end;
    end;

    @testset "Conductance" begin
        for FT in [Float32, Float64]
            vg = VanGenuchten{FT}("Loam");
            bc = BrooksCorey{FT}(vg);
            @test relative_hydraulic_conductance(vg, FT(0.2)) < 1;
            @test relative_hydraulic_conductance(bc, FT(0.2)) < 1;
            @test relative_hydraulic_conductance(vg, true, FT(-0.2)) < 1;
            @test relative_hydraulic_conductance(bc, true, FT(-0.2)) < 1;
        end;
    end;
end;
