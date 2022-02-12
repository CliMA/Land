using Photosynthesis
using Test


@testset verbose = true "Photosynthesis Test" begin
    @testset "C3 VJP" begin
        for FT in [Float32, Float64]
            leaf_1 = Leaf{FT}("C3");
            leaf_2 = Leaf{FT}("C3"; colimit = true);
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            g_mode = GCO₂Mode();
            leaf_photosynthesis!(leaf_1, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf_1, air, g_mode);
            @test true;
            leaf_photosynthesis!(leaf_2, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf_2, air, g_mode);
            @test true;
        end;
    end;

    @testset "C3 Cytochrome" begin
        for FT in [Float32, Float64]
            leaf_1 = Leaf{FT}("C3Cytochrome");
            leaf_2 = Leaf{FT}("C3Cytochrome"; colimit = true);
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            leaf_photosynthesis!(leaf_1, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf_2, air, p_mode);
            @test true;
        end;
    end;

    @testset "C4 VJP" begin
        for FT in [Float32, Float64]
            leaf_1 = Leaf{FT}("C4");
            leaf_2 = Leaf{FT}("C4"; colimit = true);
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            g_mode = GCO₂Mode();
            leaf_photosynthesis!(leaf_1, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf_1, air, g_mode);
            @test true;
            leaf_photosynthesis!(leaf_2, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf_2, air, g_mode);
            @test true;
        end;
    end;
end;
