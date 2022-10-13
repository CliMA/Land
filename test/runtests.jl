using Land
using Test


@testset verbose = true "CliMA Land" begin
    @testset "Trial" begin
        for FT in [Float32, Float32]
            vls = Land.VerticalLayers{FT}();
            Land.vertical_layers!(vls, FT(1));
            @test true;
        end
    end
end
