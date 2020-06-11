# FT and NaN test
@testset "Photosynthesis --- FT consistency and not NaN" begin
    for FT in [Float32, Float64]
        for data_set in [ PM.C3Bernacchi(FT),
                          PM.C3CLM(FT),
                          PM.C4CLM(FT),
                          PM.FluorescenceFlexas(FT),
                          PM.JmaxTDBernacchi(FT),
                          PM.JmaxTDCLM(FT),
                          PM.JmaxTDLeuning(FT),
                          PM.KcTDBernacchi(FT),
                          PM.KcTDCLM(FT),
                          PM.KoTDBernacchi(FT),
                          PM.KoTDCLM(FT),
                          PM.KpepTDBoyd(FT),
                          PM.KpepTDCLM(FT),
                          PM.RespirationTDBernacchi(FT),
                          PM.RespirationTDCLM(FT),
                          PM.VcmaxTDBernacchi(FT),
                          PM.VcmaxTDCLM(FT),
                          PM.VcmaxTDLeuning(FT),
                          PM.VomaxTDBernacchi(FT),
                          PM.VpmaxTDBoyd(FT),
                          PM.VtoRCollatz(FT),
                          PM.VtoRDefault(FT),
                          PM.ΓStarTDBernacchi(FT),
                          PM.ΓStarTDCLM(FT) ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set, FT)
        end

        mod_1     = PM.C3Bernacchi(FT)
        mod_2     = PM.C3CLM(FT)
        mod_3     = PM.C4CLM(FT)
        p_i       = FT(20.0    )
        v25       = FT(50.0    )
        j25       = FT(100.0   )
        p25       = FT(100.0   )
        par       = FT(1000.0  )
        p_O₂      = FT(21000.0 )
        t_leaf    = FT(298.3   )
        r25       = FT(1.0     )
        curvature = FT(0.9     )
        qy        = FT(0.4     )
        gsc       = FT(0.01    )
        p_atm     = FT(101325.0)
        p_a       = FT(41.0    )

        rand_T  = rand(FT) + 298
        rand_Tl = rand(FT,10) .+ 298
        for result in [ PM.arrhenius_correction(PM.KcTDBernacchi(FT), rand_T),
                        PM.arrhenius_correction(PM.JmaxTDLeuning(FT), rand_T),
                        PM.arrhenius_correction(PM.KcTDBernacchi(FT), (rand_Tl)),
                        PM.arrhenius_correction(PM.JmaxTDLeuning(FT), (rand_Tl)),
                        PM.get_jmax(PM.JmaxTDBernacchi(FT), rand(FT)*200, rand_T),
                        PM.get_kc(PM.KcTDBernacchi(FT), rand_T),
                        PM.get_ko(PM.KoTDBernacchi(FT), rand_T),
                        PM.get_ko(PM.KpepTDBoyd(FT), rand_T),
                        PM.get_r(PM.RespirationTDBernacchi(FT), rand(FT)+1, rand_T),
                        PM.get_vmax(PM.VcmaxTDBernacchi(FT), rand(FT)*100,rand_T),
                        PM.get_vmax(PM.VomaxTDBernacchi(FT), rand(FT)*100,rand_T),
                        PM.get_vmax(PM.VpmaxTDBoyd(FT), rand(FT)*100,rand_T),
                        PM.get_Γ_star(PM.ΓStarTDBernacchi(FT), rand_T),
                        PM.get_Γ_star(PM.ΓStarTDBernacchi(FT), rand_Tl),
                        PM.an_ag_r_from_pi(mod_1, p_i, v25, j25, par, p_O₂, t_leaf, r25, curvature, qy),
                        PM.an_ag_r_from_pi(mod_2, p_i, v25, j25, par, p_O₂, t_leaf, r25, curvature, qy),
                        PM.an_ag_r_from_pi(mod_3, p_i, v25, p25, par, t_leaf, r25, qy),
                        PM.an_ag_r_pi_from_gsc(mod_1, gsc, v25, j25, p25, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy),
                        PM.an_ag_r_pi_from_gsc(mod_2, gsc, v25, j25, p25, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy),
                        PM.an_ag_r_pi_from_gsc(mod_3, gsc, v25, j25, p25, p_a, t_leaf, par, p_atm, p_O₂, r25, curvature, qy),
                        PM.an_ag_r_pi_from_gsc(mod_1, gsc.*ones(FT,10), v25, j25, p25, p_a, t_leaf.*ones(FT,10), par.*ones(FT,10), p_atm, p_O₂, r25, curvature, qy),
                        PM.an_ag_r_pi_from_gsc(mod_2, gsc.*ones(FT,10), v25, j25, p25, p_a, t_leaf.*ones(FT,10), par.*ones(FT,10), p_atm, p_O₂, r25, curvature, qy),
                        PM.an_ag_r_pi_from_gsc(mod_3, gsc.*ones(FT,10), v25, j25, p25, p_a, t_leaf.*ones(FT,10), par.*ones(FT,10), p_atm, p_O₂, r25, curvature, qy),
                      ]
            recursive_FT_test(result, FT)
            recursive_NaN_test(result, FT)
        end
    end
end
