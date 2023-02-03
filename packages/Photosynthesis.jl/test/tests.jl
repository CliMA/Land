@testset verbose = true "Photosynthesis" begin
    # file temperature.jl
    @testset "Temperature" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.Arrhenius{FT}(T_REF = 298, VAL_REF = 40, ΔHA = 80000),
                        EmeraldNamespace.ArrheniusPeak{FT}(T_REF = 298, VAL_REF = 40 , ΔHA = 50000, ΔHD = 400000, ΔSV = 1000),
                        EmeraldNamespace.Q10{FT}(Q_10 = 1.4, T_REF = 298, VAL_REF = 1),
                        EmeraldNamespace.Q10Peak{FT}(Q_10 = 1.4, T_REF = 298, VAL_REF = 1, ΔHD = 400000, ΔSV = 1000)]
                @test Photosynthesis.temperature_correction(var, FT(300)) > 0;
                @test Photosynthesis.temperature_corrected_value(var, FT(300)) > 0;
                @test Photosynthesis.∂R∂T(var, FT(1), FT(300)) > 0;
            end;

            air = EmeraldNamespace.AirLayer{FT}();
            for var in [EmeraldNamespace.C3VJPModel{FT}(), EmeraldNamespace.C4VJPModel{FT}(), EmeraldNamespace.C3CytochromeModel{FT}()]
                prc = (var isa EmeraldNamespace.C3CytochromeModel ? EmeraldNamespace.CytochromeReactionCenter{FT}() : EmeraldNamespace.VJPReactionCenter{FT}();)
                Photosynthesis.photosystem_temperature_dependence!(var, prc, air, FT(300)); @test true;
                Photosynthesis.photosystem_temperature_dependence!(var, prc, air, FT(300)); @test true;
            end;

            for var in [EmeraldNamespace.Leaf{FT}(), EmeraldNamespace.Leaves1D{FT}(), EmeraldNamespace.Leaves2D{FT}()]
                @test Photosynthesis.∂R∂T(var) > 0;
            end;
        end;
    end;

    # file etr.jl, rubisco_limited.jl, light_limited.jl, product_limited.jl, and fluorescence.jl
    @testset "ETR and Rates" begin
        for FT in [Float32, Float64]
            air = EmeraldNamespace.AirLayer{FT}();
            for var in [EmeraldNamespace.Leaf{FT}(),
                        EmeraldNamespace.Leaf{FT}(PSM = EmeraldNamespace.C4VJPModel{FT}()),
                        EmeraldNamespace.Leaf{FT}(PSM = EmeraldNamespace.C3CytochromeModel{FT}(), PRC = EmeraldNamespace.CytochromeReactionCenter{FT}())]
                Photosynthesis.photosystem_temperature_dependence!(var.PSM, var.PRC, air, FT(300)); @test true;
                Photosynthesis.photosystem_electron_transport!(var.PSM, var.PRC, FT(1000), FT(20)); @test true;
                Photosynthesis.rubisco_limited_rate!(var.PSM, FT(20)); @test true;
                Photosynthesis.rubisco_limited_rate!(var.PSM, air, FT(0.1)); @test true;
                Photosynthesis.light_limited_rate!(var.PSM); @test true;
                Photosynthesis.light_limited_rate!(var.PSM, var.PRC, air, FT(0.1)); @test true;
                Photosynthesis.product_limited_rate!(var.PSM, FT(20)); @test true;
                Photosynthesis.product_limited_rate!(var.PSM, air, FT(0.1)); @test true;
                Photosynthesis.colimit_photosynthesis!(var.PSM); @test true;
                Photosynthesis.photosystem_coefficients!(var.PSM, var.PRC, FT(1000)); @test true;
            end;
        end;
    end;

    # file colimit.jl
    @testset "Colimitation" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.MinimumColimit{FT}(), EmeraldNamespace.QuadraticColimit{FT}(), EmeraldNamespace.SerialColimit{FT}(), EmeraldNamespace.SquareColimit{FT}()]
                @test Photosynthesis.colimited_rate(FT(50), FT(100), var) <= 50;
            end;

            for var in [EmeraldNamespace.C3VJPModel{FT}(), EmeraldNamespace.C4VJPModel{FT}(), EmeraldNamespace.C3CytochromeModel{FT}()]
                Photosynthesis.colimit_photosynthesis!(var); @test true;
            end;
        end;
    end;

    # file model.jl
    @testset "Core Model" begin
        for FT in [Float32, Float64]
            air = EmeraldNamespace.AirLayer{FT}();
            for var in [EmeraldNamespace.Leaf{FT}(), EmeraldNamespace.Leaves1D{FT}(), EmeraldNamespace.Leaves2D{FT}()]
                for stm in [EmeraldNamespace.AndereggSM{FT}(),
                            EmeraldNamespace.BallBerrySM{FT}(),
                            EmeraldNamespace.EllerSM{FT}(),
                            EmeraldNamespace.GentineSM{FT}(),
                            EmeraldNamespace.LeuningSM{FT}(),
                            EmeraldNamespace.MedlynSM{FT}(),
                            EmeraldNamespace.SperrySM{FT}(),
                            EmeraldNamespace.WangSM{FT}(),
                            EmeraldNamespace.Wang2SM{FT}()]
                    for bfy in [EmeraldNamespace.BetaParameterG1(), EmeraldNamespace.BetaParameterVcmax()]
                        if stm in [EmeraldNamespace.BallBerrySM{FT}(), EmeraldNamespace.LeuningSM{FT}(), EmeraldNamespace.MedlynSM{FT}()]
                            stm.β.PARAM_Y = bfy;
                        end;
                        var.SM = stm;
                        Photosynthesis.leaf_photosynthesis!(var, air, FT(0.1), FT(1000), FT(300)); @test true;
                        Photosynthesis.leaf_photosynthesis!(var, air, EmeraldNamespace.GCO₂Mode()); @test true;
                        Photosynthesis.leaf_photosynthesis!(var, air, EmeraldNamespace.PCO₂Mode()); @test true;
                    end;
                end;
            end;

            for var in [EmeraldNamespace.MonoElementSPAC{FT}(), EmeraldNamespace.MonoMLGrassSPAC{FT}(), EmeraldNamespace.MonoMLPalmSPAC{FT}(), EmeraldNamespace.MonoMLTreeSPAC{FT}()]
                Photosynthesis.leaf_photosynthesis!(var, EmeraldNamespace.GCO₂Mode()); @test true;
                Photosynthesis.leaf_photosynthesis!(var, EmeraldNamespace.PCO₂Mode()); @test true;
            end;
        end;
    end;
end;
