@testset verbose = true "EmeraldNamespace CI Coverage" begin
    # Folder air
    @testset "Air" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.AirLayer{FT}(),
                        EmeraldNamespace.Meteorology{FT}()]
                @test PkgUtility.FT_test(var, FT);
                @test PkgUtility.NaN_test(var);
            end;
        end;
    end;

    # Folder plant
    @testset "Plant" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.NonSteadyStateFlow{FT}(),
                        EmeraldNamespace.SteadyStateFlow{FT}(),
                        EmeraldNamespace.LeafHydraulics{FT}(),
                        EmeraldNamespace.RootHydraulics{FT}(),
                        EmeraldNamespace.StemHydraulics{FT}(),
                        EmeraldNamespace.BroadbandLeafBiophysics{FT}(),
                        EmeraldNamespace.HyperspectralLeafBiophysics{FT}(),
                        EmeraldNamespace.GCO₂Mode(),
                        EmeraldNamespace.PCO₂Mode(),
                        EmeraldNamespace.VDTModelAll(FT),
                        EmeraldNamespace.VDTModelDrought(FT),
                        EmeraldNamespace.VJPReactionCenter{FT}(),
                        EmeraldNamespace.CytochromeReactionCenter{FT}(),
                        EmeraldNamespace.LinearPVCurve{FT}(),
                        EmeraldNamespace.SegmentedPVCurve{FT}(),
                        EmeraldNamespace.Root{FT}(),
                        EmeraldNamespace.Stem{FT}(),
                        EmeraldNamespace.BetaParameterG1(),
                        EmeraldNamespace.BetaParameterKleaf(),
                        EmeraldNamespace.BetaParameterKsoil(),
                        EmeraldNamespace.BetaParameterPleaf(),
                        EmeraldNamespace.BetaParameterPsoil(),
                        EmeraldNamespace.BetaParameterVcmax(),
                        EmeraldNamespace.BetaParameterΘ(),
                        EmeraldNamespace.BetaFunction{FT}(),
                        EmeraldNamespace.AndereggSM{FT}(),
                        EmeraldNamespace.BallBerrySM{FT}(),
                        EmeraldNamespace.EllerSM{FT}(),
                        EmeraldNamespace.GentineSM{FT}(),
                        EmeraldNamespace.LeuningSM{FT}(),
                        EmeraldNamespace.MedlynSM{FT}(),
                        EmeraldNamespace.SperrySM{FT}(),
                        EmeraldNamespace.WangSM{FT}(),
                        EmeraldNamespace.Wang2SM{FT}(),
                        EmeraldNamespace.KcTDBernacchi(FT),
                        EmeraldNamespace.KcTDCLM(FT),
                        EmeraldNamespace.KoTDBernacchi(FT),
                        EmeraldNamespace.KoTDCLM(FT),
                        EmeraldNamespace.KpepTDCLM(FT),
                        EmeraldNamespace.KpepTDBoyd(FT),
                        EmeraldNamespace.KqTDJohnson(FT),
                        EmeraldNamespace.ΓStarTDBernacchi(FT),
                        EmeraldNamespace.ΓStarTDCLM(FT),
                        EmeraldNamespace.ηCTDJohnson(FT),
                        EmeraldNamespace.ηLTDJohnson(FT),
                        EmeraldNamespace.Q10TDAngiosperm(FT),
                        EmeraldNamespace.Q10TDGymnosperm(FT),
                        EmeraldNamespace.Q10Peak{FT}(Q_10 = 2, T_REF = 298.15, VAL_REF = 1, ΔHD = 200000, ΔSV = 500),
                        EmeraldNamespace.LogisticVC{FT}(),
                        EmeraldNamespace.PowerVC{FT}(),
                        EmeraldNamespace.WeibullVC{FT}(),
                        EmeraldNamespace.ComplexVC{FT}()]
                @test PkgUtility.FT_test(var, FT);
                @test PkgUtility.NaN_test(var);
            end;

            for var in [EmeraldNamespace.C3CytochromeModel{FT}(),
                        EmeraldNamespace.C3VJPModel{FT}(),
                        EmeraldNamespace.C4VJPModel{FT}(),
                        EmeraldNamespace.Leaf{FT}(),
                        EmeraldNamespace.Leaves1D{FT}(),
                        EmeraldNamespace.Leaves2D{FT}(),
                        EmeraldNamespace.RespirationTDBernacchi(FT),
                        EmeraldNamespace.VcmaxTDBernacchi(FT),
                        EmeraldNamespace.VomaxTDBernacchi(FT),
                        EmeraldNamespace.JmaxTDBernacchi(FT),
                        EmeraldNamespace.JmaxTDCLM(FT),
                        EmeraldNamespace.JmaxTDLeuning(FT),
                        EmeraldNamespace.RespirationTDCLM(FT),
                        EmeraldNamespace.VcmaxTDCLM(FT),
                        EmeraldNamespace.VcmaxTDLeuning(FT),
                        EmeraldNamespace.VpmaxTDBoyd(FT)]
                @test PkgUtility.FT_test(var, FT);
            end;

            # test the beta function
            tf = EmeraldNamespace.BetaFunction{FT}();
            sm = EmeraldNamespace.GentineSM{FT}();
            @test tf.FUNC(0.5) == 0.5;
            @test sm.β.FUNC(0.5) == 0.5;
        end;
    end;

    # Folder radiation
    @testset "Radiation" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.HyperspectralMLCanopyOpticalProperty{FT}(),
                        EmeraldNamespace.BroadbandSLCanopyRadiationProfile{FT}(),
                        EmeraldNamespace.HyperspectralMLCanopyRadiationProfile{FT}(),
                        EmeraldNamespace.VerhoefLIDF{FT}(),
                        EmeraldNamespace.BroadbandSLCanopy{FT}(),
                        EmeraldNamespace.HyperspectralMLCanopy{FT}(),
                        EmeraldNamespace.HyperspectralAbsorption{FT}(),
                        EmeraldNamespace.HyperspectralRadiation{FT}(),
                        EmeraldNamespace.SunSensorGeometry{FT}(),
                        EmeraldNamespace.WaveLengthSet{FT}()]
                @test PkgUtility.FT_test(var, FT);
                @test PkgUtility.NaN_test(var);
            end;
        end;
    end;

    # Folder soil
    @testset "Soil" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.BroadbandSoilAlbedo{FT}(),
                        EmeraldNamespace.HyperspectralSoilAlbedo{FT}(),
                        EmeraldNamespace.SoilLayer{FT}(),
                        EmeraldNamespace.Soil{FT}(),
                        EmeraldNamespace.BrooksCorey{FT}(500,1,"Test",1,1,1),
                        EmeraldNamespace.VanGenuchten{FT}("Sand"),
                        EmeraldNamespace.VanGenuchten{FT}("Loamy Sand"),
                        EmeraldNamespace.VanGenuchten{FT}("Sandy Loam"),
                        EmeraldNamespace.VanGenuchten{FT}("Loam"),
                        EmeraldNamespace.VanGenuchten{FT}("Sandy Clay Loam"),
                        EmeraldNamespace.VanGenuchten{FT}("Silt Loam"),
                        EmeraldNamespace.VanGenuchten{FT}("Silt"),
                        EmeraldNamespace.VanGenuchten{FT}("Clay Loam"),
                        EmeraldNamespace.VanGenuchten{FT}("Silty Clay Loam"),
                        EmeraldNamespace.VanGenuchten{FT}("Sandy Clay"),
                        EmeraldNamespace.VanGenuchten{FT}("Silty Clay"),
                        EmeraldNamespace.VanGenuchten{FT}("Clay"),
                        EmeraldNamespace.VanGenuchten{FT}("NA")]
                @test PkgUtility.FT_test(var, FT);
                @test PkgUtility.NaN_test(var);
            end;
        end;
    end;

    # Folder spac
    @testset "SPAC" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.MonoElementSPAC{FT}(),
                        EmeraldNamespace.MonoMLGrassSPAC{FT}(),
                        EmeraldNamespace.MonoMLPalmSPAC{FT}(),
                        EmeraldNamespace.MonoMLTreeSPAC{FT}()]
                @test PkgUtility.FT_test(var, FT);
            end;
        end;
    end;

    # Folder util
    @testset "Utils" begin
        for FT in [Float32, Float64]
            for var in [EmeraldNamespace.MinimumColimit{FT}(),
                        EmeraldNamespace.SerialColimit{FT}(),
                        EmeraldNamespace.SquareColimit{FT}(),
                        EmeraldNamespace.ColimitCJCLMC3(FT),
                        EmeraldNamespace.ColimitCJCLMC4(FT),
                        EmeraldNamespace.ColimitIPCLM(FT),
                        EmeraldNamespace.ColimitJCLM(FT)]
                @test PkgUtility.FT_test(var, FT);
                @test PkgUtility.NaN_test(var);
            end;
        end;
    end;
end;
