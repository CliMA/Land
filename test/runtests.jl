using ClimaCache
using PkgUtility
using Test


@testset verbose = true "ClimaCache CI Coverage" begin
    # Folder air
    @testset "Air" begin
        for FT in [Float32, Float64]
            for var in [ClimaCache.AirLayer{FT}(),
                        ClimaCache.Meteorology{FT}()]
                @test FT_test(var, FT);
                @test NaN_test(var);
            end;
        end;
    end;

    # Folder plant
    @testset "Plant" begin
        for FT in [Float32, Float64]
            for var in [ClimaCache.NonSteadyStateFlow{FT}(),
                        ClimaCache.SteadyStateFlow{FT}(),
                        ClimaCache.LeafHydraulics{FT}(),
                        ClimaCache.RootHydraulics{FT}(),
                        ClimaCache.StemHydraulics{FT}(),
                        ClimaCache.BroadbandLeafBiophysics{FT}(),
                        ClimaCache.HyperspectralLeafBiophysics{FT}(),
                        ClimaCache.GCO₂Mode(),
                        ClimaCache.PCO₂Mode(),
                        ClimaCache.VDTModelAll(FT),
                        ClimaCache.VDTModelDrought(FT),
                        ClimaCache.VJPReactionCenter{FT}(),
                        ClimaCache.CytochromeReactionCenter{FT}(),
                        ClimaCache.LinearPVCurve{FT}(),
                        ClimaCache.SegmentedPVCurve{FT}(),
                        ClimaCache.Root{FT}(),
                        ClimaCache.Stem{FT}(),
                        ClimaCache.BetaParameterG1(),
                        ClimaCache.BetaParameterKleaf(),
                        ClimaCache.BetaParameterKsoil(),
                        ClimaCache.BetaParameterPleaf(),
                        ClimaCache.BetaParameterPsoil(),
                        ClimaCache.BetaParameterVcmax(),
                        ClimaCache.BetaParameterΘ(),
                        ClimaCache.BetaFunction{FT}(),
                        ClimaCache.AndereggSM{FT}(),
                        ClimaCache.BallBerrySM{FT}(),
                        ClimaCache.EllerSM{FT}(),
                        ClimaCache.GentineSM{FT}(),
                        ClimaCache.LeuningSM{FT}(),
                        ClimaCache.MedlynSM{FT}(),
                        ClimaCache.SperrySM{FT}(),
                        ClimaCache.WangSM{FT}(),
                        ClimaCache.Wang2SM{FT}(),
                        ClimaCache.KcTDBernacchi(FT),
                        ClimaCache.KcTDCLM(FT),
                        ClimaCache.KoTDBernacchi(FT),
                        ClimaCache.KoTDCLM(FT),
                        ClimaCache.KpepTDCLM(FT),
                        ClimaCache.KpepTDBoyd(FT),
                        ClimaCache.KqTDJohnson(FT),
                        ClimaCache.ΓStarTDBernacchi(FT),
                        ClimaCache.ΓStarTDCLM(FT),
                        ClimaCache.ηCTDJohnson(FT),
                        ClimaCache.ηLTDJohnson(FT),
                        ClimaCache.Q10TDAngiosperm(FT),
                        ClimaCache.Q10TDGymnosperm(FT),
                        ClimaCache.LogisticVC{FT}(),
                        ClimaCache.PowerVC{FT}(),
                        ClimaCache.WeibullVC{FT}(),
                        ClimaCache.ComplexVC{FT}()]
                @test FT_test(var, FT);
                @test NaN_test(var);
            end;

            for var in [ClimaCache.C3CytochromeModel{FT}(),
                        ClimaCache.C3VJPModel{FT}(),
                        ClimaCache.C4VJPModel{FT}(),
                        ClimaCache.Leaf{FT}(),
                        ClimaCache.Leaves1D{FT}(),
                        ClimaCache.Leaves2D{FT}(),
                        ClimaCache.RespirationTDBernacchi(FT),
                        ClimaCache.VcmaxTDBernacchi(FT),
                        ClimaCache.VomaxTDBernacchi(FT),
                        ClimaCache.JmaxTDBernacchi(FT),
                        ClimaCache.JmaxTDCLM(FT),
                        ClimaCache.JmaxTDLeuning(FT),
                        ClimaCache.RespirationTDCLM(FT),
                        ClimaCache.VcmaxTDCLM(FT),
                        ClimaCache.VcmaxTDLeuning(FT),
                        ClimaCache.VpmaxTDBoyd(FT)]
                @test FT_test(var, FT);
            end;
        end;
    end;

    # Folder radiation
    @testset "Radiation" begin
        for FT in [Float32, Float64]
            for var in [ClimaCache.HyperspectralMLCanopyOpticalProperty{FT}(),
                        ClimaCache.BroadbandSLCanopyRadiationProfile{FT}(),
                        ClimaCache.HyperspectralMLCanopyRadiationProfile{FT}(),
                        ClimaCache.VerhoefLIDF{FT}(),
                        ClimaCache.BroadbandSLCanopy{FT}(),
                        ClimaCache.HyperspectralMLCanopy{FT}(),
                        ClimaCache.HyperspectralAbsorption{FT}(),
                        ClimaCache.HyperspectralRadiation{FT}(),
                        ClimaCache.SunSensorGeometry{FT}(),
                        ClimaCache.WaveLengthSet{FT}()]
                @test FT_test(var, FT);
                @test NaN_test(var);
            end;
        end;
    end;

    # Folder soil
    @testset "Soil" begin
        for FT in [Float32, Float64]
            for var in [ClimaCache.BroadbandSoilAlbedo{FT}(),
                        ClimaCache.HyperspectralSoilAlbedo{FT}(),
                        ClimaCache.SoilLayer{FT}(),
                        ClimaCache.Soil{FT}(),
                        ClimaCache.BrooksCorey{FT}(1,"Test",1,1,1),
                        ClimaCache.VanGenuchten{FT}("Sand"),
                        ClimaCache.VanGenuchten{FT}("Loamy Sand"),
                        ClimaCache.VanGenuchten{FT}("Sandy Loam"),
                        ClimaCache.VanGenuchten{FT}("Loam"),
                        ClimaCache.VanGenuchten{FT}("Sandy Clay Loam"),
                        ClimaCache.VanGenuchten{FT}("Silt Loam"),
                        ClimaCache.VanGenuchten{FT}("Silt"),
                        ClimaCache.VanGenuchten{FT}("Clay Loam"),
                        ClimaCache.VanGenuchten{FT}("Silty Clay Loam"),
                        ClimaCache.VanGenuchten{FT}("Sandy Clay"),
                        ClimaCache.VanGenuchten{FT}("Silty Clay"),
                        ClimaCache.VanGenuchten{FT}("Clay"),
                        ClimaCache.VanGenuchten{FT}("NA")]
                @test FT_test(var, FT);
                @test NaN_test(var);
            end;
        end;
    end;

    # Folder spac
    @testset "SPAC" begin
        for FT in [Float32, Float64]
            for var in [ClimaCache.MonoElementSPAC{FT}(),
                        ClimaCache.MonoMLGrassSPAC{FT}(),
                        ClimaCache.MonoMLPalmSPAC{FT}(),
                        ClimaCache.MonoMLTreeSPAC{FT}()]
                @test FT_test(var, FT);
            end;
        end;
    end;

    # Folder util
    @testset "Utils" begin
        constants = ClimaCache.UniversalConstants();
        @test true;

        for FT in [Float32, Float64]
            for var in [ClimaCache.MinimumColimit{FT}(),
                        ClimaCache.SerialColimit{FT}(),
                        ClimaCache.ColimitCJCLMC3(FT),
                        ClimaCache.ColimitCJCLMC4(FT),
                        ClimaCache.ColimitIPCLM(FT),
                        ClimaCache.ColimitJCLM(FT),
                        ClimaCache.AVOGADRO(FT),
                        ClimaCache.CP_D(FT),
                        ClimaCache.CP_D_MOL(FT),
                        ClimaCache.CP_I(FT),
                        ClimaCache.CP_I_MOL(FT),
                        ClimaCache.CP_L(FT),
                        ClimaCache.CP_L_MOL(FT),
                        ClimaCache.CP_V(FT),
                        ClimaCache.CP_V_MOL(FT),
                        ClimaCache.F_O₂(FT),
                        ClimaCache.GAS_R(FT),
                        ClimaCache.GRAVITY(FT),
                        ClimaCache.H_PLANCK(FT),
                        ClimaCache.K_BOLTZMANN(FT),
                        ClimaCache.K_STEFAN(FT),
                        ClimaCache.K_VON_KARMAN(FT),
                        ClimaCache.LH_V₀(FT),
                        ClimaCache.LIGHT_SPEED(FT),
                        ClimaCache.M_DRYAIR(FT),
                        ClimaCache.M_H₂O(FT),
                        ClimaCache.P_ATM(FT),
                        ClimaCache.PRESS_TRIPLE(FT),
                        ClimaCache.R_V(FT),
                        ClimaCache.RT₂₅(FT),
                        ClimaCache.T₀(FT),
                        ClimaCache.T₂₅(FT),
                        ClimaCache.T_TRIPLE(FT),
                        ClimaCache.V_H₂O(FT),
                        ClimaCache.YEAR_D(FT),
                        ClimaCache.Λ_THERMAL_H₂O(FT),
                        ClimaCache.ρ_H₂O(FT),
                        ClimaCache.ρg_MPa(FT)]
                @test FT_test(var, FT);
                @test NaN_test(var);
            end;
        end;
    end;
end;
