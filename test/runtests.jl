using ClimaCache
using StomataModels
using Test


@testset verbose = true "StomataModels Test" begin
    @testset "Conductance Limits" begin
        for FT in [Float32, Float64]
            lf_1 = ClimaCache.Leaf{FT}("C3");
            lf_2 = ClimaCache.Leaves1D{FT}("C3");
            lf_3 = ClimaCache.Leaves2D{FT}("C3");
            for lf in [lf_1, lf_2, lf_3]
                StomataModels.limit_stomatal_conductance!(lf);
                @test true;
            end;
        end;
    end;
end;
