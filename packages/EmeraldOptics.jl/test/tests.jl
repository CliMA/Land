@testset verbose = true "EmeraldOptics CI Coverage" begin
    # file photon.jl
    @testset "Photons" begin
        for FT in [Float32, Float64]
            xs = rand(FT,2);
            ys = rand(FT,2);
            EmeraldOptics.photon!(FT[400,500], xs, ys);
            @test true;
            EmeraldOptics.photon!(FT[400,500], xs);
            @test true;
            EmeraldOptics.energy!(FT[400,500], ys, xs);
            @test true;
            EmeraldOptics.energy!(FT[400,500], ys);
            @test true;
        end;
    end;

    # file transmittance.jl
    @testset "Transmittance" begin
        for FT in [Float32, Float64]
            @test EmeraldOptics.average_transmittance(FT(1), FT(1.4)) isa FT;
        end;
    end;
end;
