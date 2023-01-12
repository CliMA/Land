@testset verbose = true "StomataModels Test" begin
    @testset "Conductance Limits" begin
        for FT in [Float32, Float64]
            lf_1 = EmeraldNamespace.Leaf{FT}();
            lf_2 = EmeraldNamespace.Leaves1D{FT}();
            lf_3 = EmeraldNamespace.Leaves2D{FT}();
            for lf in [lf_1, lf_2, lf_3]
                StomataModels.limit_stomatal_conductance!(lf);
                @test true;
            end;
        end;
    end;

    @testset "∂g∂t" begin
        for FT in [Float32, Float64]
            lf_1 = EmeraldNamespace.Leaf{FT}();
            lf_2 = EmeraldNamespace.Leaves1D{FT}();
            lf_3 = EmeraldNamespace.Leaves2D{FT}();
            air  = EmeraldNamespace.AirLayer{FT}();
            for stm in [EmeraldNamespace.BallBerrySM{FT}(), EmeraldNamespace.LeuningSM{FT}(), EmeraldNamespace.MedlynSM{FT}()]
                for β_y in [EmeraldNamespace.BetaParameterG1(), EmeraldNamespace.BetaParameterVcmax()]
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
            for stm in [EmeraldNamespace.GentineSM{FT}()]
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
            for stm in [EmeraldNamespace.AndereggSM{FT}(), EmeraldNamespace.EllerSM{FT}(), EmeraldNamespace.SperrySM{FT}(), EmeraldNamespace.WangSM{FT}(), EmeraldNamespace.Wang2SM{FT}()]
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
            lf_1 = EmeraldNamespace.Leaf{FT}();
            lf_2 = EmeraldNamespace.Leaves1D{FT}();
            lf_3 = EmeraldNamespace.Leaves2D{FT}();
            air  = EmeraldNamespace.AirLayer{FT}();
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
            for spac in [EmeraldNamespace.MonoElementSPAC{FT}(),
                         EmeraldNamespace.MonoMLGrassSPAC{FT}(),
                         EmeraldNamespace.MonoMLPalmSPAC{FT}(),
                         EmeraldNamespace.MonoMLTreeSPAC{FT}()]
                StomataModels.stomatal_conductance!(spac);
                @test true;
                StomataModels.stomatal_conductance!(spac, FT(1));
                @test true;
            end;

            # TODO: add SPAC with Leaves1D in the future
            lvs = EmeraldNamespace.Leaves1D{FT}();
            air = EmeraldNamespace.AirLayer{FT}();
            StomataModels.stomatal_conductance!(lvs, air);
            @test true;
            StomataModels.stomatal_conductance!(lvs, FT(1));
            @test true;
        end;
    end;
end;
