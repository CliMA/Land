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
                kr = relative_hydraulic_conductance(vc, FT(-2.0));
                pc = critical_pressure(vc);
                @test true;
            end;
        end;
    end;

    @testset "Pressure Volume Curve" begin
        for FT in [Float32, Float64]
            pv1 = LinearPVCurve{FT}();
            pv2 = SegmentedPVCurve{FT}();
            for pv in [pv1, pv2]
                pc = xylem_pressure(pv, FT(0.8), FT(298.15));
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
                clear_legacy!(spac);
                @test true;
            end;
        end;
    end;

    @testset "Xylem End Pressure" begin
        for FT in [Float32, Float64]
            spac = MonoElementSPAC{FT}("C3");
            p = xylem_end_pressure(spac, FT(2.0));
            @test true;
            p1,p2 = xylem_end_pressure(spac, FT(1.0), FT(0.5), FT(0.5));
            @test true;
        end;
    end;
end;
