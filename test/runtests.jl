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

    @testset "PAR & APAR" begin
        for FT in [Float32, Float64]
            wls = WaveLengthSet{FT}();
            bio = LeafBiophysics{FT}(wls);
            rad = HyperspectralRadiation{FT}();

            par,apar = leaf_PAR(bio, wls, rad);
            @test true;
            par,apar = leaf_PAR(bio, wls, rad; APAR_car=false);
            @test true;
        end;
    end;

    @testset "SIF" begin
        for FT in [Float32, Float64]
            wls = WaveLengthSet{FT}();
            bio = LeafBiophysics{FT}(wls);
            rad = HyperspectralRadiation{FT}();

            sif_b,sif_f = leaf_SIF(bio, wls, rad, FT(0.01));
            @test true;
            sif_b,sif_f = leaf_SIF(bio, wls, rad, FT(0.01); ϕ_photon=false);
            @test true;
        end;
    end;

    @testset "Utils" begin
        for FT in [Float32, Float64]
            xs = rand(FT,2);
            ys = rand(FT,2);
            LeafOptics.photon!(FT[400,500], xs, ys);
            @test true;
            LeafOptics.photon!(FT[400,500], ys);
            @test true;
            LeafOptics.photon!(FT[400,500], ys, xs);
            @test true;
            LeafOptics.photon!(FT[400,500], xs);
            @test true;
        end;
    end;
end;
