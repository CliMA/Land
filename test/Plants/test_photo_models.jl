using Land

PM = Land.PhotosynthesisModels
PL = Land.Plant

@testset "PhotosynthesisModels - FT Consistency" begin
    for FT in [Float32 Float64]
        # define a list of types and then test the arrhenius correction
        for tdp in [PM.JmaxTDBernacchi{FT}(),
                    PM.JmaxTDCLM{FT}(),
                    PM.JmaxTDLeuning{FT}(),
                    PM.KcTDBernacchi{FT}(),
                    PM.KcTDCLM{FT}(),
                    PM.KoTDBernacchi{FT}(),
                    PM.KoTDCLM{FT}(),
                    PM.KpepTDBoyd{FT}(),
                    PM.KpepTDCLM{FT}(),
                    PM.RespirationTDBernacchi{FT}(),
                    PM.RespirationTDCLM{FT}(),
                    PM.VcmaxTDBernacchi{FT}(),
                    PM.VcmaxTDCLM{FT}(),
                    PM.VcmaxTDLeuning{FT}(),
                    PM.VomaxTDBernacchi{FT}(),
                    PM.VpmaxTDBoyd{FT}(),
                    PM.ΓStarTDBernacchi{FT}(),
                    PM.ΓStarTDCLM{FT}()]
            tmp = PM.arrhenius_correction(tdp, FT(298.15))
            @test tmp ≈ FT(1.0)
            @test typeof(tmp) == FT
        end

        # test the preset parameter sets, more pending
        for para_set in [PL.ESMBallBerry{FT}(),
                         PL.ESMLeuning{FT}(),
                         PL.ESMMedlyn{FT}(),
                         PM.C3VcJBernacchi{FT}(),
                         PM.C3VcVpJBernacchi{FT}(),
                         PM.C4VcVpJBoyd{FT}(),
                         PM.C4VcVpJCLM{FT}()]
            for fn in fieldnames( typeof(para_set) )
                if typeof(para_set)<:PL.AbstractEmpiricalStomatalModel
                    @test typeof( getfield(para_set,fn) ) == FT
                elseif fn != :VR
                    tmp = PM.arrhenius_correction( getfield(para_set,fn), FT(298.15) )
                    @test tmp ≈ FT(1.0)
                    @test typeof(tmp) == FT
                end
            end
        end

        # test the single return value functions
        a_net     = FT(10.0)
        curvature = FT(0.9)
        esm_1     = PL.ESMBallBerry{FT}()
        esm_2     = PL.ESMLeuning{FT}()
        esm_3     = PL.ESMMedlyn{FT}()
        g_ias_c   = FT(0.0)
        g_ias_e   = FT(0.3)
        g_list    = FT(0.1) .* ones(FT,10)
        g_max     = FT(0.8)
        gsc       = FT(0.1)
        j25       = FT(135.0)
        jtd       = PM.JmaxTDCLM{FT}()
        kctd      = PM.KcTDCLM{FT}()
        kotd      = PM.KoTDCLM{FT}()
        kpeptd    = PM.KpepTDCLM{FT}()
        mod_1     = PM.C3VcJBernacchi{FT}()
        mod_2     = PM.C3VcVpJBernacchi{FT}()
        mod_3     = PM.C4VcVpJCLM{FT}()
        p_a       = FT(40.0)
        p_atm     = FT(101325.0)
        p_H₂O     = FT(1500.0)
        p_i       = FT(10.0)
        p_O₂      = FT(21000.0)
        p25       = FT(120.0)
        par       = FT(1000.0)
        par_list  = FT(1000.0) .* ones(FT,10)
        qy        = FT(0.4)
        r25       = FT(1.5)
        r25inf    = FT(Inf)
        rh        = FT(0.5)
        rtd       = PM.RespirationTDCLM{FT}()
        t_leaf    = FT(298.15)
        t_list    = FT(298.15) .+ zeros(FT,10)
        v25       = FT(80.0)
        vpd       = FT(1500.0)
        vtd       = PM.VcmaxTDCLM{FT}()
        Γ_star    = FT(2.5)
        Γtd       = PM.ΓStarTDCLM{FT}()

        for tmp in [PM.get_j(j25, par, curvature, qy),
                    PM.get_jmax(jtd, j25, t_leaf),
                    PM.get_kc(kctd, t_leaf),
                    PM.get_ko(kotd, t_leaf),
                    PM.get_kpep(kpeptd, t_leaf),
                    PM.get_r(rtd, r25, t_leaf),
                    PM.get_r(rtd, r25inf, t_leaf),
                    PM.get_vmax(vtd, v25, t_leaf),
                    PM.get_Γ_star(Γtd, t_leaf),
                    PL.get_empirical_gsw_from_model(esm_1, a_net, p_atm, p_i, rh, vpd, Γ_star),
                    PL.get_empirical_gsw_from_model(esm_2, a_net, p_atm, p_i, rh, vpd, Γ_star),
                    PL.get_empirical_gsw_from_model(esm_3, a_net, p_atm, p_i, rh, vpd, Γ_star),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_1, esm_1),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_1, esm_2),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_1, esm_3),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_2, esm_1),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_2, esm_2),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_2, esm_3),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_3, esm_1),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_3, esm_2),
                    PL.get_empirical_gsw(curvature, g_ias_c, g_ias_e, g_max, j25, p25, p_a, p_atm, p_H₂O, p_O₂, par, qy, r25, t_leaf, v25, mod_3, esm_3)]
            @test !isnan(tmp)
            @test typeof(tmp) == FT
        end

        # test the functions for other modules
        for tmp in [PM.get_an_ag_r_from_pi(mod_1, p_i, v25, j25, par, p_O₂, t_leaf, r25, curvature, qy),
                    PM.get_an_ag_r_from_pi(mod_2, p_i, v25, j25, par, p_O₂, t_leaf, r25, curvature, qy),
                    PM.get_an_ag_r_from_pi(mod_3, p_i, v25, p25, par, t_leaf, r25, qy),
                    PM.get_an_ag_r_pi_from_gsc(mod_1, gsc, v25, j25, p25, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy),
                    PM.get_an_ag_r_pi_from_gsc(mod_2, gsc, v25, j25, p25, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy),
                    PM.get_an_ag_r_pi_from_gsc(mod_3, gsc, v25, j25, p25, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy),
                    PM.get_an_ag_r_pi_from_gsc_list(mod_1, g_list, v25, j25, p25, p_a, t_list, par_list, p_atm, p_O₂, r25, curvature, qy),
                    PM.get_an_ag_r_pi_from_gsc_list(mod_2, g_list, v25, j25, p25, p_a, t_list, par_list, p_atm, p_O₂, r25, curvature, qy),
                    PM.get_an_ag_r_pi_from_gsc_list(mod_3, g_list, v25, j25, p25, p_a, t_list, par_list, p_atm, p_O₂, r25, curvature, qy),
                    ]
            for i in tmp
                if typeof(i) <: Number
                    @test !isnan(i)
                    @test typeof(i) == FT
                else
                    for j in i
                        @test !isnan(j)
                        @test typeof(j) == FT
                    end
                end
            end
        end
    end
end