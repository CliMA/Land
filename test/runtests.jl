using ClimaCache
using PkgUtility
using Test


@testset verbose = true "ClimaCache Test" begin
    @testset "Air" begin
        for FT in [Float32, Float64]
            air = AirLayer{FT}();
            @test FT_test(air, FT);
            @test NaN_test(air);
        end;
    end;

    @testset "Plant" begin
        for FT in [Float32, Float64]
            # Leaf
            leaf_c3 = Leaf{FT}("C3");
            leaf_c4 = Leaf{FT}("C4");
            leaf_cy = Leaf{FT}("C3Cytochrome");
            wls = WaveLengthSet{FT}(collect(400:10:2500));
            leaf_d3 = Leaf{FT}("C3", wls);
            leaf_d4 = Leaf{FT}("C4", wls);
            leaf_dy = Leaf{FT}("C3Cytochrome", wls);
            for leaf in [leaf_c3, leaf_c4, leaf_cy, leaf_d3, leaf_d4, leaf_dy]
                @test FT_test(leaf, FT);
                # NaN test will not pass because of the NaNs in temperature dependency structures
                # @test NaN_test(leaf);
            end;

            # LeafBiophysics
            lbio1 = LeafBiophysics{FT}();
            lbio2 = LeafBiophysics{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
            for lbio in [lbio1, lbio2]
                @test FT_test(lbio, FT);
                @test NaN_test(lbio);
            end;

            # Fluorescence model
            vdt1 = VanDerTolFluorescenceModel{FT}();
            vdt2 = VanDerTolFluorescenceModel{FT}(true);
            cyto = CytochromeFluorescenceModel{FT}();
            for vdt in [vdt1, vdt2, cyto]
                @test FT_test(vdt, FT);
                @test NaN_test(vdt);
            end;

            # Reaction center
            rc1 = VJPReactionCenter{FT}();
            rc2 = CytochromeReactionCenter{FT}();
            for rc in [rc1, rc2]
                @test FT_test(rc, FT);
                @test NaN_test(rc);
            end;

            # Photosynthesis model
            cy_1 = C3CytochromeModel{FT}();
            cy_2 = C3CytochromeModel{FT}(v_cmax25 = 30, r_d25 = 1);
            c3_1 = C3VJPModel{FT}();
            c3_2 = C3VJPModel{FT}(v_cmax25 = 30, j_max25 = 50, r_d25 = 1);
            c4_1 = C4VJPModel{FT}();
            c4_2 = C4VJPModel{FT}(v_cmax25 = 30, v_pmax25 = 40, r_d25 = 1);
            for st in [cy_1, cy_2, c3_1, c3_2, c4_1, c4_2]
                for rc in [rc1, rc2]
                    @test FT_test(st, FT);
                    # NaN test will not pass because of the NaNs in temperature dependency structures
                    # @test NaN_test(st);
                end;
            end;

            # Mode and colimitations
            mod1 = GCO₂Mode();
            mod2 = PCO₂Mode();
            col1 = MinimumColimit{FT}();
            col2 = QuadraticColimit{FT}(0.98);
            for st in [mod1, mod2, col1, col2]
                for rc in [rc1, rc2]
                    @test FT_test(st, FT);
                    @test NaN_test(st);
                end;
            end;

            # Temperature dependency
            td_1 = Arrhenius{FT}(298.15, 41.0, 79430.0);
            td_2 = ArrheniusPeak{FT}(298.15, 1.0, 57500.0, 439000.0, 1400.0);
            td_3 = Q10{FT}(298.15, 0.0140/8760, 1.4);
            for td in [td_1, td_2, td_3]
                @test FT_test(td, FT);
                @test NaN_test(td);
            end;
        end;
    end;

    @testset "Radiation" begin
        for FT in [Float32, Float64]
            wls1 = WaveLengthSet{FT}();
            wls2 = WaveLengthSet{FT}(collect(400:5:2500));
            wls3 = WaveLengthSet{FT}(collect(400:5:2500); opti=ClimaCache.OPTI_2017);
            for wls in [wls1, wls2, wls3]
                @test FT_test(wls, FT);
                # NaN test will not pass because of the NaNs in wls2 and wls3
                # @test NaN_test(wls);
            end;
        end;
    end;
end;
