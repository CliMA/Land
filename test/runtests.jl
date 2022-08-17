using ClimaCache
using LeafOptics
using Test


@testset verbose = true "LeafOptics CI Coverage" begin
    # file spectra.jl
    @testset "Spectra" begin
        for FT in [Float32, Float64]
            wls  = ClimaCache.WaveLengthSet{FT}();
            bio  = ClimaCache.HyperspectralLeafBiophysics{FT}();
            lha  = ClimaCache.HyperspectralAbsorption{FT}();
            spac = ClimaCache.MonoMLTreeSPAC{FT}();

            leaf_spectra!(bio, wls, lha, FT(50));
            @test true;
            leaf_spectra!(bio, wls, lha, FT(50));
            @test true;
            leaf_spectra!(bio, wls, lha, FT(49); APAR_car = false);
            @test true;
            leaf_spectra!(bio, wls, lha, FT(48); reabsorb = false);
            @test true;
            leaf_spectra!(bio, wls, FT(0.1), FT(0.45), FT(0.05), FT(0.25));
            @test true;
            leaf_spectra!(spac);
            @test true;
        end;
    end;

    # file radiation.jl
    @testset "PAR & APAR" begin
        for FT in [Float32, Float64]
            wls = ClimaCache.WaveLengthSet{FT}();
            bio = ClimaCache.HyperspectralLeafBiophysics{FT}();
            rad = ClimaCache.HyperspectralRadiation{FT}();

            par,apar,ppar = leaf_PAR(bio, wls, rad);
            @test true;
            par,apar,ppar = leaf_PAR(bio, wls, rad; APAR_car=false);
            @test true;
        end;
    end;

    # file fluorescence.jl
    @testset "SIF" begin
        for FT in [Float32, Float64]
            wls = ClimaCache.WaveLengthSet{FT}();
            bio = ClimaCache.HyperspectralLeafBiophysics{FT}();
            rad = ClimaCache.HyperspectralRadiation{FT}();

            sif_b,sif_f = leaf_SIF(bio, wls, rad, FT(0.01));
            @test true;
            sif_b,sif_f = leaf_SIF(bio, wls, rad, FT(0.01); ϕ_photon = false);
            @test true;
        end;
    end;

    # file photon.jl
    @testset "Utils" begin
        for FT in [Float32, Float64]
            xs = rand(FT,2);
            ys = rand(FT,2);
            LeafOptics.photon!(FT[400,500], xs, ys);
            @test true;
            LeafOptics.photon!(FT[400,500], xs);
            @test true;
            LeafOptics.energy!(FT[400,500], ys, xs);
            @test true;
            LeafOptics.energy!(FT[400,500], ys);
            @test true;
        end;
    end;
end;
