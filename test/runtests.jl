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

    @testset "1D+2D Leaves" begin
        for FT in [Float32, Float64]
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            g_mode = GCO₂Mode();
            for LT in ["C3", "C4", "C3Cytochrome"]
                leaves_1d = Leaves1D{FT}(LT);
                leaves_2d = Leaves2D{FT}(LT);
                leaf_photosynthesis!(leaves_1d, air, g_mode);
                @test true
                leaf_photosynthesis!(leaves_1d, air, p_mode);
                @test true
                leaf_photosynthesis!(leaves_2d, air, g_mode);
                @test true
                leaf_photosynthesis!(leaves_2d, air, p_mode);
                @test true
            end;
        end;
    end

    @testset "P&G Modes" begin
        for FT in [Float32, Float64]
            air    = AirLayer{FT}();
            p_mode = PCO₂Mode();
            g_mode = GCO₂Mode();
            for LT in ["C3", "C4", "C3Cytochrome"]
                leaf_1 = Leaf{FT}(LT);
                leaf_2 = Leaf{FT}(LT);
                for glc in collect(0.05:0.05:0.3)
                    leaf_photosynthesis!(leaf_1, air, g_mode, FT(glc));
                    leaf_photosynthesis!(leaf_2, air, p_mode, leaf_1.p_CO₂_i);
                    @test leaf_1.PSM.a_gross ≈ leaf_2.PSM.a_gross;
                    @test leaf_1.PSM.e_to_c ≈ leaf_2.PSM.e_to_c;
                    @test leaf_1.PRC.ϕ_f ≈ leaf_2.PRC.ϕ_f;
                end;
            end;
        end;
    end;
end;
