@testset verbose = true "EmeraldConstants Test" begin
    @testset "Structures" begin
        EmeraldConstants.UNIVERSAL_CONSTANTS.CP_D = 1000;
        @test EmeraldConstants.UNIVERSAL_CONSTANTS.CP_D == 1000;
    end;

    @testset "Wrapper functions" begin
        for FT in [Float32, Float64]
            for var in [EmeraldConstants.AVOGADRO(FT),
                        EmeraldConstants.CP_D(FT),
                        EmeraldConstants.CP_D_MOL(FT),
                        EmeraldConstants.CP_I(FT),
                        EmeraldConstants.CP_I_MOL(FT),
                        EmeraldConstants.CP_L(FT),
                        EmeraldConstants.CP_L_MOL(FT),
                        EmeraldConstants.CP_V(FT),
                        EmeraldConstants.CP_V_MOL(FT),
                        EmeraldConstants.F_O₂(FT),
                        EmeraldConstants.GAS_R(FT),
                        EmeraldConstants.GRAVITY(FT),
                        EmeraldConstants.H_PLANCK(FT),
                        EmeraldConstants.K_BOLTZMANN(FT),
                        EmeraldConstants.K_STEFAN(FT),
                        EmeraldConstants.K_VON_KARMAN(FT),
                        EmeraldConstants.LH_V₀(FT),
                        EmeraldConstants.LIGHT_SPEED(FT),
                        EmeraldConstants.M_DRYAIR(FT),
                        EmeraldConstants.M_H₂O(FT),
                        EmeraldConstants.P_ATM(FT),
                        EmeraldConstants.PRESS_TRIPLE(FT),
                        EmeraldConstants.R_V(FT),
                        EmeraldConstants.RT₂₅(FT),
                        EmeraldConstants.T₀(FT),
                        EmeraldConstants.T₂₅(FT),
                        EmeraldConstants.T_TRIPLE(FT),
                        EmeraldConstants.V_H₂O(FT),
                        EmeraldConstants.YEAR_D(FT),
                        EmeraldConstants.Λ_THERMAL_H₂O(FT),
                        EmeraldConstants.ρ_H₂O(FT),
                        EmeraldConstants.ρg_MPa(FT)]
                @test !isnan(var);
                @test typeof(var) == FT;
            end;
        end;
    end;
end;
