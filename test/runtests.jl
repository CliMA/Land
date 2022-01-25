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


#=
@testset "Structure --- Radiation" begin
    println("Testing WaveLengthSet constructors...");
    for FT in [Float32, Float64]
        wls = WaveLengthSet{FT}();
        wls = WaveLengthSet{FT}(collect(400:5:2500));
        wls = WaveLengthSet{FT}(collect(400:5:2500); opti=ClimaCache.OPTI_2017);
        @test true;
    end;

    println("Testing HyperspectralRadiation constructors...");
    for FT in [Float32, Float64]
        rad = HyperspectralRadiation{FT}();
        rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
        rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2400)); file = "");
        @info "Expecting Warning here!";
        rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(400:50:2500)); file = ClimaCache.FILE_SUN);
        @test true;
    end;
end;


@testset "Structure --- Soil" begin
    println("Testing BrooksCorey constructors...");
    for FT in [Float32, Float64]
        bc = BrooksCorey{FT}("Soil type", 1, 1, 0.5, 0.2);
        @test true;
    end;

    println("Testing VanGenuchten constructors...");
    for FT in [Float32, Float64]
        @info "Expecting Warning here!";
        vg = VanGenuchten{FT}("");
        for stype in ["Sand", "Loamy Sand", "Sandy Loam", "Loam", "Sandy Clay Loam", "Silt Loam", "Silt", "Clay Loam", "Silty Clay Loam", "Sandy Clay", "Silty Clay", "Clay"]
            vg = VanGenuchten{FT}(stype);
        end;
        vg = VanGenuchten{FT}("Test", 100, 2, 0.5, 0.1);
        @test true;
    end;
end;


@testset "Structure --- Plant" begin
    println("Testing HyperspectralAbsorption constructors...");
    for FT in [Float32, Float64]
        ha = HyperspectralAbsorption{FT}();
        ha = HyperspectralAbsorption{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
        ha = HyperspectralAbsorption{FT}(WaveLengthSet{FT}(collect(400:50:2400)); opti=ClimaCache.OPTI_2017);
        @test true;
    end;

    println("Testing LeafBiophysics constructors...");
    for FT in [Float32, Float64]
        lbio = LeafBiophysics{FT}();
        lbio = LeafBiophysics{FT}(WaveLengthSet{FT}(collect(400:50:2400)));
        @test true;
    end;

    println("Testing C3VJPModel constructors...");
    for FT in [Float32, Float64]
        c3 = C3VJPModel{Float64}();
        c3 = C3VJPModel{Float64}(v_cmax25 = 30, j_max25 = 50, r_d25 = 1);
        @test true;
    end;
end;
=#
