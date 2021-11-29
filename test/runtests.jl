using ClimaCache
using Test


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
        vg = VanGenuchten{FT}("Loam");
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

    println("Testing C₃VJPSystem constructors...");
    for FT in [Float32, Float64]
        c3 = C₃VJPSystem{Float64}();
        c3 = C₃VJPSystem{Float64}(v_max25 = 30, j_max25 = 50, r_d25 = 1);
        @test true;
    end;
end;
