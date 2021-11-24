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
        ha = HyperspectralAbsorption{FT}(WaveLengthSet{FT}(collect(400:50:2500)); opti=ClimaCache.OPTI_2017);
        @test true;
    end;
end;
