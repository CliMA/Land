using ClimaCache
using PlantHydraulics
using Test


@testset verbose = true "PlantHydraulics Test" begin
    @testset "Vulnerability Curve" begin
        for FT in [Float32, Float64]
            vc1 = LogisticVC{FT}(2, 2);
            vc2 = PowerVC{FT}(2, 2);
            vc3 = WeibullVC{FT}(2, 2);
            vs1 = ComplexVC{FT}([0.3,0.3,0.4], [vc1, vc1, vc1]);
            vs2 = ComplexVC{FT}([0.3,0.3,0.4], [vc2, vc2, vc2]);
            vs3 = ComplexVC{FT}([0.3,0.3,0.4], [vc3, vc3, vc3]);
            vs4 = ComplexVC{FT}([0.3,0.3,0.4], [vc1, vc2, vc3]);
            for vc in [vc1, vc2, vc3, vs1, vs2, vs3, vs4]
                kr = PlantHydraulics.relative_hydraulic_conductance(vc, FT(-2.0));
                pc = PlantHydraulics.critical_pressure(vc);
                @test true;
            end;
        end;
    end;

    @testset "Pressure Volume Curve" begin
        for FT in [Float32, Float64]
            pv1 = LinearPVCurve{FT}();
            pv2 = SegmentedPVCurve{FT}();
            for pv in [pv1, pv2]
                pc = PlantHydraulics.xylem_pressure(pv, FT(0.8), FT(298.15));
                @test true;
            end;
        end;
    end;

    @testset "Cavitation Legacy" begin
        for FT in [Float32, Float64]
            spac1 = MonoElementSPAC{FT}("C3");
            spac2 = MonoGrassSPAC{FT}("C3");
            spac3 = MonoPalmSPAC{FT}("C3");
            spac4 = MonoTreeSPAC{FT}("C3");
            for spac in [spac1, spac2, spac3, spac4]
                PlantHydraulics.clear_legacy!(spac);
                @test true;
            end;
        end;
    end;

    @testset "Xylem End Pressure" begin
        for FT in [Float32, Float64]
            spac = MonoElementSPAC{FT}("C3");
            p = PlantHydraulics.xylem_end_pressure(spac, FT(2.0));
            @test true;
            p1,p2 = PlantHydraulics.xylem_end_pressure(spac, FT(1.0), FT(0.5), FT(0.5));
            @test true;
        end;
    end;

    @testset "Flow Profile" begin
        for FT in [Float32, Float64]
            spac1 = MonoElementSPAC{FT}("C3");
            spac2 = MonoGrassSPAC{FT}("C3");
            spac3 = MonoPalmSPAC{FT}("C3");
            spac4 = MonoTreeSPAC{FT}("C3");
            spac5 = MonoElementSPAC{FT}("C3"; ssm = false);
            spac6 = MonoGrassSPAC{FT}("C3"; ssm = false);
            spac7 = MonoPalmSPAC{FT}("C3"; ssm = false);
            spac8 = MonoTreeSPAC{FT}("C3"; ssm = false);
            spacx = MonoTreeSPAC{FT}("C3"; ssm = false);
            spacy = MonoTreeSPAC{FT}("C3"; ssm = false);
            spacx.ROOTS[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacx.BRANCHES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacx.LEAVES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.ROOTS[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.BRANCHES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.LEAVES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.TRUNK.HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            for spac in [spac1, spac2, spac3, spac4, spac5, spac6, spac7, spac8, spacx, spacy]
                xylem_flow_profile!(spac, FT(10));
                @test true;
            end;
        end;
    end;

    @testset "Pressure Profile" begin
        for FT in [Float32, Float64]
            spac1 = MonoElementSPAC{FT}("C3");
            spac2 = MonoGrassSPAC{FT}("C3");
            spac3 = MonoPalmSPAC{FT}("C3");
            spac4 = MonoTreeSPAC{FT}("C3");
            spac5 = MonoElementSPAC{FT}("C3"; ssm = false);
            spac6 = MonoGrassSPAC{FT}("C3"; ssm = false);
            spac7 = MonoPalmSPAC{FT}("C3"; ssm = false);
            spac8 = MonoTreeSPAC{FT}("C3"; ssm = false);
            spacx = MonoTreeSPAC{FT}("C3"; ssm = false);
            spacy = MonoTreeSPAC{FT}("C3"; ssm = false);
            spacx.ROOTS[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacx.BRANCHES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacx.LEAVES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.ROOTS[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.BRANCHES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.LEAVES[1].HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            spacy.TRUNK.HS.FLOW = ClimaCache.SteadyStateFlow{FT}(0);
            for spac in [spac1, spac2, spac3, spac4, spac5, spac6, spac7, spac8, spacx, spacy]
                xylem_pressure_profile!(spac);
                @test true;
            end;
        end;
    end;

    @testset "Critical Flow" begin
        for FT in [Float32, Float64]
            lhs = LeafHydraulics{FT}();
            @test PlantHydraulics.critical_flow(lhs, FT(298.15)) > 0;
            spac = MonoElementSPAC{FT}("C3");
            @test PlantHydraulics.critical_flow(spac) > 0;
        end;
    end;
end;
