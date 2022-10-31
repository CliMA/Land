using ClimaCache
using PlantHydraulics
using Test


@testset verbose = true "PlantHydraulics Test" begin
    @testset "Vulnerability Curve" begin
        for FT in [Float32, Float64]
            vc1 = ClimaCache.LogisticVC{FT}();
            vc2 = ClimaCache.PowerVC{FT}();
            vc3 = ClimaCache.WeibullVC{FT}();
            vs1 = ClimaCache.ComplexVC{FT}([0.3,0.3,0.4], [vc1, vc1, vc1]);
            vs2 = ClimaCache.ComplexVC{FT}([0.3,0.3,0.4], [vc2, vc2, vc2]);
            vs3 = ClimaCache.ComplexVC{FT}([0.3,0.3,0.4], [vc3, vc3, vc3]);
            vs4 = ClimaCache.ComplexVC{FT}([0.3,0.3,0.4], [vc1, vc2, vc3]);
            for vc in [vc1, vc2, vc3, vs1, vs2, vs3, vs4]
                kr = PlantHydraulics.relative_hydraulic_conductance(vc, FT(-2.0));
                pc = PlantHydraulics.critical_pressure(vc);
                @test true;
            end;
        end;
    end;

    @testset "Pressure Volume Curve" begin
        for FT in [Float32, Float64]
            pv1 = ClimaCache.LinearPVCurve{FT}();
            pv2 = ClimaCache.SegmentedPVCurve{FT}();
            for pv in [pv1, pv2]
                pc = PlantHydraulics.xylem_pressure(pv, FT(0.8), FT(298.15));
                @test true;
            end;
        end;
    end;

    @testset "Cavitation Legacy" begin
        for FT in [Float32, Float64]
            spac1 = ClimaCache.MonoElementSPAC{FT}();
            spac2 = ClimaCache.MonoMLGrassSPAC{FT}();
            spac3 = ClimaCache.MonoMLPalmSPAC{FT}();
            spac4 = ClimaCache.MonoMLTreeSPAC{FT}();
            for spac in [spac1, spac2, spac3, spac4]
                PlantHydraulics.clear_legacy!(spac);
                @test true;
            end;
        end;
    end;

    @testset "Xylem End Pressure" begin
        for FT in [Float32, Float64]
            spac = ClimaCache.MonoElementSPAC{FT}();
            p = PlantHydraulics.xylem_end_pressure(spac, FT(2.0));
            @test true;
            # p1,p2 = PlantHydraulics.xylem_end_pressure(spac, FT(1.0), FT(0.5), FT(0.5));
            # @test true;
        end;
    end;

    @testset "Flow Profile" begin
        for FT in [Float32, Float64]
            spac1 = ClimaCache.MonoElementSPAC{FT}();
            spac2 = ClimaCache.MonoMLGrassSPAC{FT}();
            spac3 = ClimaCache.MonoMLPalmSPAC{FT}();
            spac4 = ClimaCache.MonoMLTreeSPAC{FT}();
            spacx = ClimaCache.MonoMLTreeSPAC{FT}();
            spacy = ClimaCache.MonoMLTreeSPAC{FT}();
            spacx.ROOTS[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacx.BRANCHES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacx.LEAVES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 1);
            spacy.ROOTS[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacy.BRANCHES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacy.LEAVES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 1);
            spacy.TRUNK.HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            for spac in [spac1, spac2, spac3, spac4, spacx, spacy]
                PlantHydraulics.xylem_flow_profile!(spac, FT(10));
                @test true;
            end;
        end;
    end;

    @testset "Pressure Profile" begin
        for FT in [Float32, Float64]
            spac1 = ClimaCache.MonoElementSPAC{FT}();
            spac2 = ClimaCache.MonoMLGrassSPAC{FT}();
            spac3 = ClimaCache.MonoMLPalmSPAC{FT}();
            spac4 = ClimaCache.MonoMLTreeSPAC{FT}();
            spacx = ClimaCache.MonoMLTreeSPAC{FT}();
            spacy = ClimaCache.MonoMLTreeSPAC{FT}();
            spacx.ROOTS[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacx.BRANCHES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacx.LEAVES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 1);
            spacy.ROOTS[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacy.BRANCHES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            spacy.LEAVES[1].HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 1);
            spacy.TRUNK.HS.FLOW = ClimaCache.NonSteadyStateFlow{FT}(DIM_CAPACITY = 5);
            for spac in [spac1, spac2, spac3, spac4, spacx, spacy]
                PlantHydraulics.xylem_pressure_profile!(spac);
                @test true;
            end;
        end;
    end;

    @testset "∂E∂P" begin
        for FT in [Float32, Float64]
            leaf1 = ClimaCache.Leaf{FT}();
            leaf2 = ClimaCache.Leaves1D{FT}();
            leaf3 = ClimaCache.Leaves2D{FT}();
            PlantHydraulics.∂E∂P(leaf1, FT(0));
            @test true;
            PlantHydraulics.∂E∂P(leaf2, FT(0), 1);
            @test true;
            PlantHydraulics.∂E∂P(leaf3, FT(0));
            @test true;
        end;
    end;

    @testset "Critical Flow" begin
        for FT in [Float32, Float64]
            lhs = ClimaCache.LeafHydraulics{FT}();
            @test PlantHydraulics.critical_flow(lhs, FT(298.15)) > 0;
            spac = ClimaCache.MonoElementSPAC{FT}();
            @test PlantHydraulics.critical_flow(spac) > 0;
        end;
    end;

    @testset "Tuning Factor" begin
        for FT in [Float32, Float64]
            spac1 = ClimaCache.MonoElementSPAC{FT}();
            spac2 = ClimaCache.MonoMLGrassSPAC{FT}();
            spac3 = ClimaCache.MonoMLPalmSPAC{FT}();
            spac4 = ClimaCache.MonoMLTreeSPAC{FT}();
            for spac in [spac1, spac2, spac3, spac4]
                PlantHydraulics.β_factor!(spac);
                @test true;
            end;
        end;
    end;

    @testset "Energy Budget" begin
        for FT in [Float32, Float64]
            spac1 = ClimaCache.MonoMLGrassSPAC{FT}();
            spac2 = ClimaCache.MonoMLPalmSPAC{FT}();
            spac3 = ClimaCache.MonoMLTreeSPAC{FT}();
            for spac in [spac1, spac2, spac3]
                PlantHydraulics.plant_energy!(spac);
                @test true;
                PlantHydraulics.plant_energy!(spac, FT(1));
                @test true;
            end;
        end;
    end;
end;
