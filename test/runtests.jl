using ClimaCache
using StomataModels
using Test


@testset verbose = true "StomataModels Test" begin
    @testset "Conductance Limits" begin
        for FT in [Float32, Float64]
            lf_1 = ClimaCache.Leaf{FT}();
            lf_2 = ClimaCache.Leaves1D{FT}();
            lf_3 = ClimaCache.Leaves2D{FT}();
            for lf in [lf_1, lf_2, lf_3]
                StomataModels.limit_stomatal_conductance!(lf);
                @test true;
            end;
        end;
    end;

    @testset "∂g∂t" begin
        for FT in [Float32, Float64]
            lf_1 = ClimaCache.Leaf{FT}();
            lf_2 = ClimaCache.Leaves1D{FT}();
            lf_3 = ClimaCache.Leaves2D{FT}();
            air  = ClimaCache.AirLayer{FT}();
            for stm in [ClimaCache.BallBerrySM{FT}(), ClimaCache.LeuningSM{FT}(), ClimaCache.MedlynSM{FT}()]
                for β_y in [ClimaCache.BetaParameterG1(), ClimaCache.BetaParameterVcmax()]
                    stm.β.PARAM_Y = β_y;
                    lf_1.SM = stm;
                    lf_2.SM = stm;
                    lf_3.SM = stm;
                    StomataModels.∂g∂t(lf_1, air);
                    @test true;
                    StomataModels.∂g∂t(lf_2, air, 1);
                    @test true;
                    StomataModels.∂g∂t(lf_3, air);
                    StomataModels.∂g∂t(lf_3, air, 1);
                    @test true;
                end;
            end;
            for stm in [ClimaCache.GentineSM{FT}()]
                lf_1.SM = stm;
                lf_2.SM = stm;
                lf_3.SM = stm;
                StomataModels.∂g∂t(lf_1, air);
                @test true;
                StomataModels.∂g∂t(lf_2, air, 1);
                @test true;
                StomataModels.∂g∂t(lf_3, air);
                StomataModels.∂g∂t(lf_3, air, 1);
                @test true;
            end;
            for stm in [ClimaCache.AndereggSM{FT}(), ClimaCache.EllerSM{FT}(), ClimaCache.SperrySM{FT}(), ClimaCache.WangSM{FT}(), ClimaCache.Wang2SM{FT}()]
                lf_1.SM = stm;
                lf_2.SM = stm;
                lf_3.SM = stm;
                StomataModels.∂g∂t(lf_1, air);
                @test true;
                StomataModels.∂g∂t(lf_2, air, 1);
                @test true;
                StomataModels.∂g∂t(lf_3, air);
                StomataModels.∂g∂t(lf_3, air, 1);
                @test true;
            end;
        end;
    end;

    @testset "∂gₙ∂t" begin
        for FT in [Float32, Float64]
            lf_1 = ClimaCache.Leaf{FT}();
            lf_2 = ClimaCache.Leaves1D{FT}();
            lf_3 = ClimaCache.Leaves2D{FT}();
            air  = ClimaCache.AirLayer{FT}();
            StomataModels.∂gₙ∂t(lf_1, air);
            @test true;
            StomataModels.∂gₙ∂t(lf_2, air);
            @test true;
            StomataModels.∂gₙ∂t(lf_3, air);
            @test true;
        end;
    end;

    @testset "Prognostic Conductance" begin
        for FT in [Float32, Float64]
            for spac in [ClimaCache.MonoElementSPAC{FT}(),
                         ClimaCache.MonoMLGrassSPAC{FT}(),
                         ClimaCache.MonoMLPalmSPAC{FT}(),
                         ClimaCache.MonoMLTreeSPAC{FT}()]
                StomataModels.stomatal_conductance!(spac);
                @test true;
                StomataModels.stomatal_conductance!(spac, FT(1));
                @test true;
            end;

            # TODO: add SPAC with Leaves1D in the future
            lvs = ClimaCache.Leaves1D{FT}();
            air = ClimaCache.AirLayer{FT}();
            StomataModels.stomatal_conductance!(lvs, air);
            @test true;
            StomataModels.stomatal_conductance!(lvs, air, FT(1));
            @test true;
        end;
    end;
end;
