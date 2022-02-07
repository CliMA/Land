using Photosynthesis
using Test


@testset verbose = true "Photosynthesis Test" begin
    @testset "C3 VJP" begin
        for FT in [Float32, Float64]
            leaf   = Leaf{FT}("C3");
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            g_mode = GCO₂Mode();
            leaf_photosynthesis!(leaf, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf, air, g_mode);
            @test true;
        end;
    end;

    @testset "C3 Cytochrome" begin
        for FT in [Float32, Float64]
            leaf   = Leaf{FT}("C3Cytochrome");
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            leaf_photosynthesis!(leaf, air, p_mode);
            @test true;
        end;
    end;

    @testset "C4 VJP" begin
        for FT in [Float32, Float64]
            leaf   = Leaf{FT}("C4");
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            g_mode = GCO₂Mode();
            leaf_photosynthesis!(leaf, air, p_mode);
            @test true;
            leaf_photosynthesis!(leaf, air, g_mode);
            @test true;
        end;
    end;
end;
