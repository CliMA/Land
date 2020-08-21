# Function tests adapted from experiments/Radiation_test_BRDF
println("\nTesting the layered model...");
@testset "CanopyLayers --- RT function test" begin
    for FT in [Float32, Float64]
        canopy_rt = create_canopy_rt(FT, nLayer=20);
        wl_set    = create_wave_length(FT);
        rt_dim    = create_rt_dims(canopy_rt, wl_set);
        canRad_rt = create_canopy_rads(FT, rt_dim);
        canOpt_rt = create_canopy_opticals(FT, rt_dim);
        sunRad_rt = create_incoming_radiation(wl_set);
        soil      = create_soil_opticals(wl_set);
        angles    = SolarAngles{FT}();
        rt_con    = create_rt_container(FT, rt_dim);

        leaf_1    = create_leaf_bios(FT, rt_dim);
        leaf_2    = create_leaf_bios(FT, rt_dim);

        collections = initialize_rt_module(FT, LAI=FT(3));

        wl_blue   = FT(450.0);
        wl_red    = FT(600.0);
        wl_FarRed = FT(740.0);
        wl_Red    = FT(685.0);
        ind_wle_blue = argmin( abs.(wl_set.WLE .- wl_blue  ) );
        ind_wle_red  = argmin( abs.(wl_set.WLE .- wl_red   ) );
        ind_wlf_FR   = argmin( abs.(wl_set.WLF .- wl_FarRed) );
        ind_wlf_R    = argmin( abs.(wl_set.WLF .- wl_Red   ) );
        ind_red      = argmin( abs.(wl_set.WL  .- wl_Red   ) );
        ind_NIR      = argmin( abs.(wl_set.WL  .- 800      ) );

        leaf_2.Cab = FT(80.0 );
        leaf_2.Cw  = FT(0.012);
        fluspect!(leaf_1, wl_set);
        fluspect!(leaf_2, wl_set);

        arrayOfLeaves = [create_leaf_bios(FT, rt_dim) for i in 1:rt_dim.nLayer];
        for i in 1:rt_dim.nLayer
            fluspect!(arrayOfLeaves[i], wl_set);
        end

        canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
        canopy_matrices!(arrayOfLeaves, canOpt_rt);
        short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
        canopy_fluxes!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, arrayOfLeaves, wl_set, rt_con);

        SIF_FR  = FT[];
        SIF_R   = FT[];
        reflVIS = FT[];
        reflNIR = FT[];

        angles.tts    = 30;
        angles.psi    = 0;
        canopy_rt.LAI = 3;
        VZA           = collect(FT, -89.5:0.5:89.5);

        for VZA_ in VZA
            angles.tto = VZA_;
            canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
            canopy_matrices!(arrayOfLeaves, canOpt_rt);
            short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
            SIF_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set, rt_con, rt_dim);

            push!(reflVIS, canRad_rt.alb_obs[ind_red   ]);
            push!(reflNIR, canRad_rt.alb_obs[ind_NIR   ]);
            push!(SIF_R  , canRad_rt.SIF_obs[ind_wlf_R ]);
            push!(SIF_FR , canRad_rt.SIF_obs[ind_wlf_FR]);
        end

        reflVIS = FT[];
        reflNIR = FT[];
        SIF_FR  = FT[];
        SIF_R   = FT[];

        angles.tts    = 48;
        angles.psi    = 0;
        canopy_rt.LAI = FT(3.22);

        for psi in 0:5:360
            angles.psi = psi;
            for VZA in 0:5:85
                angles.tto = VZA;

                canopy_geometry!(canopy_rt, angles, canOpt_rt, rt_con);
                canopy_matrices!(arrayOfLeaves, canOpt_rt);
                short_wave!(canopy_rt, canOpt_rt, canRad_rt, sunRad_rt, soil, rt_con);
                SIF_fluxes!(arrayOfLeaves, canOpt_rt, canRad_rt, canopy_rt, soil, wl_set, rt_con, rt_dim);

                push!(reflVIS, canRad_rt.alb_obs[28]);
                push!(reflNIR, canRad_rt.alb_obs[52]);
                push!(SIF_R  , canRad_rt.SIF_obs[8 ]);
                push!(SIF_FR , canRad_rt.SIF_obs[20]);
            end
        end

        A         = reshape(reflNIR, (18,73));
        B         = reshape(reflVIS, (18,73));
        SIFFER    = reshape(SIF_R  , (18,73));
        SIFFER_FR = reshape(SIF_FR , (18,73));

        # do something?
        @test true;
    end
end
