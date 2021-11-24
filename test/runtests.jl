using ClimaCache
using Test

@testset "Structure --- Radiation" begin
    println("             Testing WaveLengthSet constructors...");
    for FT in [Float32, Float64]
        wls = WaveLengthSet{FT}();
        wls = WaveLengthSet{FT}(collect(FT,400:5:2500));
        wls = WaveLengthSet{FT}(collect(FT,400:5:2500); opti=ClimaCache.OPTI_2017);
        @test true;
    end;

    println("             Testing HyperspectralRadiation constructors...");
    for FT in [Float32, Float64]
        rad = HyperspectralRadiation{FT}();
        rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(FT,400:50:2400)));
        rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(FT,400:50:2400)), "");
        rad = HyperspectralRadiation{FT}(WaveLengthSet{FT}(collect(FT,400:50:2400)), ClimaCache.FILE_SUN);
        @test true;
    end;
end
