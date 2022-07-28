using ClimaCache
using Photosynthesis
using Test


@testset verbose = true "Photosynthesis CI Coverage" begin
    # file temperature.jl
    @testset "Temperature" begin
        for FT in [Float32, Float64]
            for var in [ClimaCache.Arrhenius{FT}(T_REF = 298, VAL_REF = 40, ΔHA = 80000),
                        ClimaCache.ArrheniusPeak{FT}(T_REF = 298, VAL_REF = 40 , ΔHA = 50000, ΔHD = 400000, ΔSV = 1000),
                        ClimaCache.Q10{FT}(Q_10 = 1.4, T_REF = 298, VAL_REF = 1)]
                @test Photosynthesis.temperature_correction(var, FT(300)) > 0;
                @test Photosynthesis.temperature_corrected_value(var, FT(300)) > 0;
                @test Photosynthesis.∂R∂T(var, FT(1), FT(300)) > 0;
            end;

            air = ClimaCache.AirLayer{FT}();
            for var in [ClimaCache.C3VJPModel{FT}(), ClimaCache.C4VJPModel{FT}(), ClimaCache.C3CytochromeModel{FT}()]
                Photosynthesis.photosystem_temperature_dependence!(var, air, FT(300)); @test true;
                Photosynthesis.photosystem_temperature_dependence!(var, air, FT(300)); @test true;
            end;

            for var in [ClimaCache.Leaf{FT}(), ClimaCache.Leaves1D{FT}(), ClimaCache.Leaves2D{FT}()]
                @test Photosynthesis.∂R∂T(var) > 0;
            end;
        end;
    end;

    # file etr.jl, rubisco_limited.jl, light_limited.jl, product_limited.jl, and fluorescence.jl
    @testset "ETR and Rates" begin
        for FT in [Float32, Float64]
            air = ClimaCache.AirLayer{FT}();
            for var in [ClimaCache.Leaf{FT}(),
                        ClimaCache.Leaf{FT}(PSM = ClimaCache.C4VJPModel{FT}()),
                        ClimaCache.Leaf{FT}(PSM = ClimaCache.C3CytochromeModel{FT}(), PRC = ClimaCache.CytochromeReactionCenter{FT}())]
                Photosynthesis.photosystem_temperature_dependence!(var.PSM, air, FT(300)); @test true;
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
            for var in [ClimaCache.MinimumColimit{FT}(), ClimaCache.QuadraticColimit{FT}(), ClimaCache.SerialColimit{FT}()]
                @test Photosynthesis.colimited_rate(FT(50), FT(100), var) <= 50;
            end;

            for var in [ClimaCache.C3VJPModel{FT}(), ClimaCache.C4VJPModel{FT}(), ClimaCache.C3CytochromeModel{FT}()]
                Photosynthesis.colimit_photosynthesis!(var); @test true;
            end;
        end;
    end;

    # file model.jl
    @testset "Core Model" begin
        for FT in [Float32, Float64]
            air = ClimaCache.AirLayer{FT}();
            for var in [ClimaCache.Leaf{FT}(), ClimaCache.Leaves1D{FT}(), ClimaCache.Leaves2D{FT}()]
                for stm in [ClimaCache.AndereggSM{FT}(),
                            ClimaCache.BallBerrySM{FT}(),
                            ClimaCache.EllerSM{FT}(),
                            ClimaCache.GentineSM{FT}(),
                            ClimaCache.LeuningSM{FT}(),
                            ClimaCache.MedlynSM{FT}(),
                            ClimaCache.SperrySM{FT}(),
                            ClimaCache.WangSM{FT}(),
                            ClimaCache.Wang2SM{FT}()]
                    for bfy in [ClimaCache.BetaParameterG1(), ClimaCache.BetaParameterVcmax()]
                        if stm in [ClimaCache.BallBerrySM{FT}(), ClimaCache.LeuningSM{FT}(), ClimaCache.MedlynSM{FT}()]
                            stm.β.PARAM_Y = bfy;
                        end;
                        var.SM = stm;
                        leaf_photosynthesis!(var, air, FT(0.1), FT(1000), FT(300)); @test true;
                        leaf_photosynthesis!(var, air, ClimaCache.GCO₂Mode()); @test true;
                        leaf_photosynthesis!(var, air, ClimaCache.PCO₂Mode()); @test true;
                    end;
                end;
            end;
        end;
    end;
end;
