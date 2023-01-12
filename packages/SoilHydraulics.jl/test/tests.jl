@testset verbose = true "SoilHydraulics Test" begin
    @testset "Constructor" begin
        for FT in [Float32, Float64]
            bc = EmeraldNamespace.BrooksCorey{FT}(EmeraldNamespace.VanGenuchten{FT}("Loam"));
            @test true;
        end;
    end;

    @testset "Potential" begin
        for FT in [Float32, Float64]
            vg = EmeraldNamespace.VanGenuchten{FT}("Loam");
            bc = EmeraldNamespace.BrooksCorey{FT}(vg);
            @test SoilHydraulics.soil_ψ_25(vg, FT(0.2)) < 0;
            @test SoilHydraulics.soil_ψ_25(bc, FT(0.2)) < 0;
            @test SoilHydraulics.soil_ψ_25(vg, FT(1.0)) == 0;
            @test SoilHydraulics.soil_ψ_25(bc, FT(1.0)) == -bc.Ψ_SAT;
        end;
    end;

    @testset "Content" begin
        for FT in [Float32, Float64]
            vg = EmeraldNamespace.VanGenuchten{FT}("Loam");
            bc = EmeraldNamespace.BrooksCorey{FT}(vg);
            @test vg.Θ_RES < SoilHydraulics.soil_θ(vg, FT(-1)) < vg.Θ_SAT;
            @test SoilHydraulics.soil_θ(vg, FT(0)) == vg.Θ_SAT;
            @test bc.Θ_RES < SoilHydraulics.soil_θ(bc, FT(-1)) < bc.Θ_SAT;
            @test SoilHydraulics.soil_θ(bc, FT(0)) == bc.Θ_SAT;
        end;
    end;

    @testset "Conductance" begin
        for FT in [Float32, Float64]
            vg = EmeraldNamespace.VanGenuchten{FT}("Loam");
            bc = EmeraldNamespace.BrooksCorey{FT}(vg);
            @test SoilHydraulics.relative_hydraulic_conductance(vg, FT(0.2)) < 1;
            @test SoilHydraulics.relative_hydraulic_conductance(bc, FT(0.2)) < 1;
            @test SoilHydraulics.relative_hydraulic_conductance(vg, true, FT(-0.2)) < 1;
            @test SoilHydraulics.relative_hydraulic_conductance(bc, true, FT(-0.2)) < 1;
            @test SoilHydraulics.relative_hydraulic_conductance(vg, true, FT(0.01)) == 1;
            @test SoilHydraulics.relative_hydraulic_conductance(bc, true, FT(0.01)) == 1;
        end;
    end;

    @testset "Budgets" begin
        for FT in [Float32, Float64]
            spac = EmeraldNamespace.MonoMLTreeSPAC{FT}();
            SoilHydraulics.soil_budget!(spac);
            @test true;
            SoilHydraulics.soil_budget!(spac, FT(1));
            @test true;
        end;
    end;
end;
