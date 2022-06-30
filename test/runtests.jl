using ClimaCache
using PkgUtility
using Test


@testset verbose = true "ClimaCache Test" begin
    @testset "Soil" begin
        for FT in [Float32, Float64]
            soil1 = ClimaCache.Soil{FT}(FT[0,-1]);
            soil2 = ClimaCache.Soil{FT}(FT[0,-1], true);
            for soil in [soil1, soil2]
                @test FT_test(soil, FT);
                @test NaN_test(soil);
            end;
        end;
    end;

    @testset "Air" begin
        for FT in [Float32, Float64]
            air = ClimaCache.AirLayer{FT}();
            @test FT_test(air, FT);
            @test NaN_test(air);
        end;
    end;

    @testset "Plant" begin
        for FT in [Float32, Float64]
            # Pressure volume curve
            pvc1 = ClimaCache.LinearPVCurve{FT}();
            pvc2 = ClimaCache.SegmentedPVCurve{FT}();
            for pvc in [pvc1, pvc2]
                @test FT_test(pvc, FT);
                @test NaN_test(pvc);
            end;

            # Plant hydraulic system
            lhs1 = ClimaCache.LeafHydraulics{FT}();
            lhs2 = ClimaCache.LeafHydraulics{FT}(ssm = true);
            lhs3 = ClimaCache.LeafHydraulics{FT}(ssm = false);
            rhs1 = ClimaCache.RootHydraulics{FT}();
            rhs2 = ClimaCache.RootHydraulics{FT}(ssm = true);
            rhs3 = ClimaCache.RootHydraulics{FT}(ssm = false);
            shs1 = ClimaCache.StemHydraulics{FT}();
            shs2 = ClimaCache.StemHydraulics{FT}(ssm = true);
            shs3 = ClimaCache.StemHydraulics{FT}(ssm = false);
            for hs in [lhs1, lhs2, lhs3, rhs1, rhs2, rhs3, shs1, shs2, shs3]
                @test FT_test(hs, FT);
                @test NaN_test(hs);
            end;

            # Root and Stem
            root = ClimaCache.Root{FT}();
            stem = ClimaCache.Stem{FT}();
            for hs in [root, stem]
                @test FT_test(hs, FT);
                @test NaN_test(hs);
            end;

            # Leaf
            leaf_c3 = ClimaCache.Leaf{FT}("C3");
            leaf_c4 = ClimaCache.Leaf{FT}("C4");
            leaf_cy = ClimaCache.Leaf{FT}("C3Cytochrome");
            leaf_d3 = ClimaCache.Leaf{FT}("C3", colimit = true);
            leaf_d4 = ClimaCache.Leaf{FT}("C4", colimit = true);
            leaf_dy = ClimaCache.Leaf{FT}("C3Cytochrome", colimit = true);
            wls     = ClimaCache.WaveLengthSet{FT}(collect(400:10:2500));
            leaf_e3 = ClimaCache.Leaf{FT}("C3", wls);
            leaf_e4 = ClimaCache.Leaf{FT}("C4", wls);
            leaf_ey = ClimaCache.Leaf{FT}("C3Cytochrome", wls);
            for leaf in [leaf_c3, leaf_c4, leaf_cy, leaf_d3, leaf_d4, leaf_dy, leaf_e3, leaf_e4, leaf_ey]
                @test FT_test(leaf, FT);
                # NaN test will not pass because of the NaNs in temperature dependency structures
                # @test NaN_test(leaf);
            end;

            # Leaves1D
            leaves_c3 = ClimaCache.Leaves1D{FT}("C3");
            leaves_c4 = ClimaCache.Leaves1D{FT}("C4");
            leaves_cy = ClimaCache.Leaves1D{FT}("C3Cytochrome");
            leaves_d3 = ClimaCache.Leaves1D{FT}("C3", colimit = true);
            leaves_d4 = ClimaCache.Leaves1D{FT}("C4", colimit = true);
            leaves_dy = ClimaCache.Leaves1D{FT}("C3Cytochrome", colimit = true);
            for leaves in [leaves_c3, leaves_c4, leaves_cy, leaves_d3, leaves_d4, leaves_dy]
                @test FT_test(leaves, FT);
                # NaN test will not pass because of the NaNs in temperature dependency structures
                # @test NaN_test(leaf);
            end;

            # Leaves2D
            leaves_c3 = ClimaCache.Leaves2D{FT}("C3");
            leaves_c4 = ClimaCache.Leaves2D{FT}("C4");
            leaves_cy = ClimaCache.Leaves2D{FT}("C3Cytochrome");
            leaves_d3 = ClimaCache.Leaves2D{FT}("C3", colimit = true);
            leaves_d4 = ClimaCache.Leaves2D{FT}("C4", colimit = true);
            leaves_dy = ClimaCache.Leaves2D{FT}("C3Cytochrome", colimit = true);
            wls       = ClimaCache.WaveLengthSet{FT}(collect(400:10:2500));
            leaves_e3 = ClimaCache.Leaves2D{FT}("C3", wls);
            leaves_e4 = ClimaCache.Leaves2D{FT}("C4", wls);
            leaves_ey = ClimaCache.Leaves2D{FT}("C3Cytochrome", wls);
            for leaves in [leaves_c3, leaves_c4, leaves_cy, leaves_d3, leaves_d4, leaves_dy, leaves_e3, leaves_e4, leaves_ey]
                @test FT_test(leaves, FT);
                # NaN test will not pass because of the NaNs in temperature dependency structures
                # @test NaN_test(leaf);
            end;

            # LeafBiophysics
            lbio1 = ClimaCache.HyperspectralLeafBiophysics{FT}();
            lbio2 = ClimaCache.HyperspectralLeafBiophysics{FT}(ClimaCache.WaveLengthSet{FT}(collect(400:50:2400)));
            lbio3 = ClimaCache.BroadbandLeafBiophysics{FT}();
            for lbio in [lbio1, lbio2, lbio3]
                @test FT_test(lbio, FT);
                @test NaN_test(lbio);
            end;

            # Fluorescence model
            vdt1 = ClimaCache.VanDerTolFluorescenceModel{FT}();
            vdt2 = ClimaCache.VanDerTolFluorescenceModel{FT}(true);
            for vdt in [vdt1, vdt2]
                @test FT_test(vdt, FT);
                @test NaN_test(vdt);
            end;

            # Reaction center
            rc1 = ClimaCache.VJPReactionCenter{FT}();
            rc2 = ClimaCache.CytochromeReactionCenter{FT}();
            for rc in [rc1, rc2]
                @test FT_test(rc, FT);
                @test NaN_test(rc);
            end;

            # Photosynthesis model
            cy_1 = ClimaCache.C3CytochromeModel{FT}();
            cy_2 = ClimaCache.C3CytochromeModel{FT}(v_cmax25 = 30, r_d25 = 1, colimit = true);
            c3_1 = ClimaCache.C3VJPModel{FT}();
            c3_2 = ClimaCache.C3VJPModel{FT}(v_cmax25 = 30, j_max25 = 50, r_d25 = 1, colimit = true);
            c4_1 = ClimaCache.C4VJPModel{FT}();
            c4_2 = ClimaCache.C4VJPModel{FT}(v_cmax25 = 30, v_pmax25 = 40, r_d25 = 1, colimit = true);
            for st in [cy_1, cy_2, c3_1, c3_2, c4_1, c4_2]
                for rc in [rc1, rc2]
                    @test FT_test(st, FT);
                    # NaN test will not pass because of the NaNs in temperature dependency structures
                    # @test NaN_test(st);
                end;
            end;

            # Mode and colimitations
            mod1 = ClimaCache.GCO₂Mode();
            mod2 = ClimaCache.PCO₂Mode();
            col1 = ClimaCache.MinimumColimit{FT}();
            col2 = ClimaCache.QuadraticColimit{FT}(0.98);
            col3 = ClimaCache.SerialColimit{FT}();
            for st in [mod1, mod2, col1, col2, col3]
                for rc in [rc1, rc2]
                    @test FT_test(st, FT);
                    @test NaN_test(st);
                end;
            end;

            # Temperature dependency
            td_1 = ClimaCache.Arrhenius{FT}(298.15, 41.0, 79430.0);
            td_2 = ClimaCache.ArrheniusPeak{FT}(298.15, 1.0, 57500.0, 439000.0, 1400.0);
            td_3 = ClimaCache.Q10{FT}(298.15, 0.0140/8760, 1.4);
            for td in [td_1, td_2, td_3]
                @test FT_test(td, FT);
                @test NaN_test(td);
            end;

            # Beta function
            β_1 = ClimaCache.BetaFunctionLeafK(ClimaCache.BetaParameterG1());
            β_2 = ClimaCache.BetaFunctionLeafK(ClimaCache.BetaParameterVcmax());
            β_3 = ClimaCache.BetaFunctionSoilK(ClimaCache.BetaParameterG1());
            β_4 = ClimaCache.BetaFunctionSoilK(ClimaCache.BetaParameterVcmax());
            for td in [β_1, β_2, β_3, β_4]
                @test FT_test(td, FT);
                @test NaN_test(td);
            end;
        end;
    end;

    @testset "Radiation" begin
        for FT in [Float32, Float64]
            # Wave length sets
            wls1 = ClimaCache.WaveLengthSet{FT}();
            wls2 = ClimaCache.WaveLengthSet{FT}(collect(400:5:2500));
            wls3 = ClimaCache.WaveLengthSet{FT}(collect(400:5:2500); opti=ClimaCache.OPTI_2017);
            for wls in [wls1, wls2, wls3]
                @test FT_test(wls, FT);
                # NaN test will not pass because of the NaNs in wls2 and wls3
                # @test NaN_test(wls);
            end;

            # Solar radiation
            rad1 = ClimaCache.BroadbandRadiation{FT}();
            rad2 = ClimaCache.HyperspectralRadiation{FT}();
            for rad in [rad1, rad2]
                @test FT_test(rad, FT);
                @test NaN_test(rad);
            end;

            # Sun-sensor geometry
            ssg = ClimaCache.SunSensorGeometry{FT}();
            @test FT_test(ssg, FT);
            @test NaN_test(ssg);

            # Canopy structure
            can = ClimaCache.HyperspectralMLCanopy{FT}();
            @test FT_test(can, FT);
            @test NaN_test(can);
        end;
    end;

    @testset "SPAC" begin
        for FT in [Float32, Float64]
            spac1 = ClimaCache.MonoElementSPAC{FT}("C3");
            spac2 = ClimaCache.MonoMLGrassSPAC{FT}("C3");
            spac3 = ClimaCache.MonoMLPalmSPAC{FT}("C3");
            spac4 = ClimaCache.MonoMLTreeSPAC{FT}("C3");
            spac5 = ClimaCache.MonoElementSPAC{FT}("C3"; ssm = false);
            spac6 = ClimaCache.MonoMLGrassSPAC{FT}("C3"; ssm = false);
            spac7 = ClimaCache.MonoMLPalmSPAC{FT}("C3"; ssm = false);
            spac8 = ClimaCache.MonoMLTreeSPAC{FT}("C3"; ssm = false);
            for spac in [spac1, spac2, spac3, spac4, spac5, spac6, spac7, spac8]
                @test FT_test(spac, FT);
                # NaN test will not pass because of the NaNs in temperature dependency structures
                # @test NaN_test(wls);
            end;
        end;
    end;
end;
