# FT and NaN tests
# TODO test the MArrays as well
@testset "CanopyRT --- FT and NaN test" begin
    for FT in [Float32, Float64]
        for data_set in [ RT.Canopy4RT{FT, 20, 3.0}(),
                          RT.create_canopy_optical(FT, 10, 10, 10, 10; using_marray=false),
                          #RT.create_canopy_optical(FT, 10, 10, 10, 10; using_marray=true),
                          RT.CanopyRadiation{FT, 10, 10, 10, 10, 20}(),
                          RT.create_incoming_radiation(FT, FT[600,620,640]; using_marray=false),
                          #RT.create_incoming_radiation(FT, FT[600,620,640]; using_marray=true),
                          RT.create_leaf_bio(FT, 10, 10, 10; using_marray=false),
                          #RT.create_leaf_bio(FT, 10, 10, 10; using_marray=true),
                          RT.create_opti_par(FT, FT[600,620,640]; using_marray=false),
                          #RT.create_opti_par(FT, FT[600,620,640]; using_marray=true),
                          RT.SoilOpti{FT}(FT[600,610,620], FT[0.1,0.1,0.1], FT[0.1,0.1,0.1], FT(280.0)),
                          RT.SolarAngles{FT}(),
                          RT.create_wl_para_set(FT)
                           ]
            recursive_FT_test(data_set, FT)
            recursive_NaN_test(data_set, FT)
        end
    end
end




# Function tests adapted from experiments/Radiation_test_BRDF
@testset "CanopyRT --- RT function test" begin
    for FT in [Float32, Float64]
        wl_set    = RT.create_wl_para_set(FT)
        leaf_1    = RT.create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF)
        leaf_2    = RT.create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF)
        canopy_rt = RT.Canopy4RT{FT, 20, 3.0}()
        canRad_rt = RT.CanopyRadiation{FT, wl_set.nwl, wl_set.nWlF, length(canopy_rt.litab), length(canopy_rt.lazitab), canopy_rt.nlayers}()
        canOpt_rt = RT.create_canopy_optical(FT, wl_set.nwl, canopy_rt.nlayers, length(canopy_rt.lazitab), length(canopy_rt.litab); using_marray=false)
        sunRad_rt = RT.create_incoming_radiation(FT, wl_set.swl)
        soil      = RT.SoilOpti{FT}(wl_set.wl, FT(0.2)*ones(FT, length(wl_set.wl)), FT[0.1], FT(290.0))
        angles    = RT.SolarAngles{FT}()

        wl_blue   = FT(450.0)
        wl_red    = FT(600.0)
        wl_FarRed = FT(740.0)
        wl_Red    = FT(685.0)
        ind_wle_blue = argmin( abs.(wl_set.wle .- wl_blue  ) )
        ind_wle_red  = argmin( abs.(wl_set.wle .- wl_red   ) )
        ind_wlf_FR   = argmin( abs.(wl_set.wlf .- wl_FarRed) )
        ind_wlf_R    = argmin( abs.(wl_set.wlf .- wl_Red   ) )
        ind_red      = argmin( abs.(wl_set.wl  .- wl_Red   ) )
        ind_NIR      = argmin( abs.(wl_set.wl  .- 800      ) )

        leaf_2.Cab = FT(80.0 )
        leaf_2.Cw  = FT(0.012)
        RT.fluspect!(leaf_1, canopy_rt, wl_set)
        RT.fluspect!(leaf_2, canopy_rt, wl_set)

        arrayOfLeaves = [RT.create_leaf_bio(FT, wl_set.nwl, wl_set.nWlE, wl_set.nWlF) for i in 1:canopy_rt.nlayers]
        for i in 1:canopy_rt.nlayers
            RT.fluspect!(arrayOfLeaves[i], canopy_rt, wl_set)
        end

        RT.compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
        RT.compute_canopy_matrices!(arrayOfLeaves, canOpt_rt)
        RT.simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil)
        RT.derive_canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set)

        SIF_FR  = FT[]
        SIF_R   = FT[]
        reflVIS = FT[]
        reflNIR = FT[]

        angles.tts    = 30
        angles.psi    = 0
        canopy_rt.LAI = 3
        VZA           = collect(FT, -89.5:0.5:89.5)

        for VZA_ in VZA
            angles.tto = VZA_
            RT.compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
            RT.compute_canopy_matrices!(arrayOfLeaves, canOpt_rt)
            RT.simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil)
            RT.computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set)
            
            push!(reflVIS, canRad_rt.alb_obs[ind_red   ])
            push!(reflNIR, canRad_rt.alb_obs[ind_NIR   ])
            push!(SIF_R  , canRad_rt.SIF_obs[ind_wlf_R ])
            push!(SIF_FR , canRad_rt.SIF_obs[ind_wlf_FR])
        end

        reflVIS = FT[]
        reflNIR = FT[]
        SIF_FR  = FT[]
        SIF_R   = FT[]

        angles.tts    = 48
        angles.psi    = 0
        canopy_rt.LAI = FT(3.22)

        for psi in 0:360
            angles.psi = psi
            for VZA in 0:1:85
                angles.tto = VZA

                RT.compute_canopy_geometry!(canopy_rt, angles, canOpt_rt)
                RT.compute_canopy_matrices!(arrayOfLeaves, canOpt_rt)
                RT.simulate_short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil)
                RT.computeSIF_Fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set)

                push!(reflVIS, canRad_rt.alb_obs[28])
                push!(reflNIR, canRad_rt.alb_obs[52])
                push!(SIF_R  , canRad_rt.SIF_obs[8 ])
                push!(SIF_FR , canRad_rt.SIF_obs[20])
            end
        end

        A         = reshape(reflNIR, (86,361))
        B         = reshape(reflVIS, (86,361))
        SIFFER    = reshape(SIF_R  , (86,361))
        SIFFER_FR = reshape(SIF_FR , (86,361))

        # do something?
        
        @test true
    end
end
