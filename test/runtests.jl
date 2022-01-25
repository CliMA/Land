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
