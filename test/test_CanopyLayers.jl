# Test big leaf canopy model
println("\nTesting the big leaf model...");
@testset "CanopyLayers --- big leaf model" begin
    for FT in [Float32, Float64]
        for result in [ big_leaf_partition(FT(3.0), FT(30.0), FT(1000.0)),
                        big_leaf_partition(FT(2.0), FT(30.0), FT(1000.0)),
                        big_leaf_partition(FT(1.0), FT(30.0), FT(1000.0)) ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end




# FT and NaN tests
println("\nTesting the SCOPE model...");
@testset "CanopyLayers --- SCOPE model" begin
    for FT in [Float32, Float64]
        collections = initialize_rt_module(FT; nLayer=20, LAI=3);
        angles, can, can_opt, can_rad, in_rad, leaves, rt_con, rt_dim, soil, wls = collections;
        for data_set in collections
            @test FT_test(data_set, FT);
            @test NaN_test(data_set);
        end

        # add more tests
        angles.tts = 30;
        angles.psi = 0;
        can.LAI    = 3;
        VZA        = collect(FT, -89.5:0.5:89.5);
        for VZA_ in VZA
            angles.tto = VZA_;
            canopy_geometry!(can, angles, can_opt, rt_con);
            canopy_matrices!(leaves, can_opt);
            short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
            SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);
        end
        @test true;

        angles.tts = 48;
        angles.psi = 0;
        can.LAI    = FT(3.22);
        for psi in 0:5:360
            angles.psi = psi;
            for VZA in 0:5:85
                angles.tto = VZA;
                canopy_geometry!(can, angles, can_opt, rt_con);
                canopy_matrices!(leaves, can_opt);
                short_wave!(can, can_opt, can_rad, in_rad, soil, rt_con);
                SIF_fluxes!(leaves, can_opt, can_rad, can, soil, wls, rt_con, rt_dim);
            end
        end
        @test true;

        # test warnings
        println("Expect warnings here!");
        warn_wls   = create_wave_length(FT, collect(FT,2100:100:2600));
        warn_inrad = create_incoming_radiation(warn_wls);
        @test true;

        # add more tests that has not been used
        can.clump_a = 0.6;
        can.clump_b = 0.1;
        CanopyLayers.clumping_factor!(can, angles);
        canopy_fluxes!(can, can_opt, can_rad, in_rad, soil, leaves[1:1], wls, rt_con);
        canopy_matrices!(leaves[1:1], can_opt);
        SIF_fluxes!(leaves[1:1], can_opt, can_rad, can, soil, wls, rt_con, rt_dim);
        thermal_fluxes!(leaves[1:1], can_opt, can_rad, can, soil, [FT(400.0)], wls);
        thermal_fluxes!(leaves[1:2], can_opt, can_rad, can, soil, [FT(400.0)], wls);
        @test true;

        # utility functions
        CanopyLayers.dcum(FT(1.1), FT(1.1), FT(1.1));
        CanopyLayers.e2phot(rand(FT,10), rand(FT,10));
        CanopyLayers.volscatt!(rand(FT,4), FT(40), FT(90), FT(40), FT(0));
        @test true;
    end
end




# indices test
println("\nTesting the Indicies...");
@testset "CanopyLayers --- Indicies" begin
    for FT in [Float32, Float64]
        collections = initialize_rt_module(FT; nLayer=20, LAI=3);
        angles, can, can_opt, can_rad, in_rad, leaves, rt_con, rt_dim, soil, wls = collections;
        _indices = [EVI(can_rad, wls),
                    EVI2(can_rad, wls),
                    LSWI(can_rad, wls),
                    NDVI(can_rad, wls),
                    NIRv(can_rad, wls),
                    SIF_740(can_rad, wls),
                    SIF_757(can_rad, wls),
                    SIF_771(can_rad, wls)];
        @test NaN_test(_indices);
        @test FT_test(_indices, FT);

        # out bounds warnings
        println("Expect warnings here!");
        REF_WL(can_rad, wls, FT(3000));
        SIF_WL(can_rad, wls, FT(3000));
        @test true
    end
end
