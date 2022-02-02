using LeafOptics
using Test


@testset verbose = true "LeafOptics Test" begin
    @testset "Spectra" begin
        for FT in [Float32, Float64]
            wls = WaveLengthSet{FT}();
            bio = LeafBiophysics{FT}(wls);
            lha = HyperspectralAbsorption{FT}(wls);

            leaf_spectra!(bio, wls, lha);
            @test true;
            leaf_spectra!(bio, wls, lha; APAR_car=false);
            @test true;
            leaf_spectra!(bio, wls, lha; APAR_car=false, α=FT(59));
            @test true;

            leaf_spectra!(bio, wls, FT(0.1), FT(0.45), FT(0.05), FT(0.25));
            @test true;
        end;
    end;
end;
