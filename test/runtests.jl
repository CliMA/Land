using ClimaCache
using Photosynthesis
using Test


@testset verbose = true "Photosynthesis Test" begin
    @testset "C3 VJP" begin
        for FT in [Float32, Float64]
            leaf   = ClimaCache.Leaf{FT}();
            air    = ClimaCache.AirLayer{FT}();
            p_mode = ClimaCache.PCO₂Mode();
            g_mode = ClimaCache.GCO₂Mode();
            leaf_photosynthesis!(leaf, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf, air, g_mode);
            @test true;
        end;
    end;

    @testset "C3 Cytochrome" begin
        for FT in [Float32, Float64]
            leaf   = ClimaCache.Leaf{FT}(PSM = ClimaCache.C3CytochromeModel{FT}(), PRC = ClimaCache.CytochromeReactionCenter{FT}());
            air    = ClimaCache.AirLayer{FT}();
            p_mode = ClimaCache.PCO₂Mode();
            g_mode = ClimaCache.GCO₂Mode();
            leaf_photosynthesis!(leaf, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf, air, g_mode);
            @test true;
        end;
    end;

    @testset "C4 VJP" begin
        for FT in [Float32, Float64]
            leaf   = ClimaCache.Leaf{FT}(PSM = ClimaCache.C4VJPModel{FT}());
            air    = ClimaCache.AirLayer{FT}();
            p_mode = ClimaCache.PCO₂Mode();
            g_mode = ClimaCache.GCO₂Mode();
            leaf_photosynthesis!(leaf, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf, air, g_mode);
            @test true;
        end;
    end;

    @testset "1D+2D Leaves" begin
        for FT in [Float32, Float64]
            air    = ClimaCache.AirLayer{FT}();
            p_mode = ClimaCache.PCO₂Mode();
            g_mode = ClimaCache.GCO₂Mode();
            leaves_1d = ClimaCache.Leaves1D{FT}();
            leaves_2d = ClimaCache.Leaves2D{FT}();
            leaf_photosynthesis!(leaves_1d, air, g_mode);
            @test true;
            leaf_photosynthesis!(leaves_1d, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaves_2d, air, g_mode);
            @test true;
            leaf_photosynthesis!(leaves_2d, air, p_mode);
            @test true;
        end;
    end

    @testset "P&G Modes" begin
        for FT in [Float32, Float64]
            air    = ClimaCache.AirLayer{FT}();
            p_mode = ClimaCache.PCO₂Mode();
            g_mode = ClimaCache.GCO₂Mode();
            leaf_1 = ClimaCache.Leaf{FT}();
            leaf_2 = ClimaCache.Leaf{FT}();
            for glc in collect(0.05:0.05:0.3)
                leaf_1._g_CO₂ = glc;
                leaf_photosynthesis!(leaf_1, air, g_mode);
                leaf_2._p_CO₂_i = leaf_1._p_CO₂_i;
                leaf_photosynthesis!(leaf_2, air, p_mode);
                @test leaf_1.PSM.a_gross ≈ leaf_2.PSM.a_gross;
                @test leaf_1.PSM._e_to_c ≈ leaf_2.PSM._e_to_c;
                @test leaf_1.PRC.ϕ_f ≈ leaf_2.PRC.ϕ_f;
            end;
        end;
    end;

    @testset "SPAC" begin
        for FT in [Float32, Float64]
            p_mode = ClimaCache.PCO₂Mode();
            g_mode = ClimaCache.GCO₂Mode();
            spac1 = ClimaCache.MonoElementSPAC{FT}();
            spac2 = ClimaCache.MonoMLGrassSPAC{FT}();
            spac3 = ClimaCache.MonoMLPalmSPAC{FT}();
            spac4 = ClimaCache.MonoMLTreeSPAC{FT}();
            for spac in [spac1, spac2, spac3, spac4]
                leaf_photosynthesis!(spac, g_mode);
                @test true;
                leaf_photosynthesis!(spac, p_mode);
                @test true;
            end;
        end;
    end;
end;
