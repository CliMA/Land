@testset verbose = true "PlantHydraulics" begin
    @testset "Structures" begin
        for FT in [Float32, Float64]
            leaf  = LeafHydraulics{FT}();
            root  = RootHydraulics{FT}();
            stem  = StemHydraulics{FT}();
            grass = create_grass(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
            palm  = create_palm(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            tree  = create_tree(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            treet = TreeSimple{FT}();
            _vc1  = WeibullSingle{FT}();
            _vc2  = WeibullDual{FT}();
            _sh1  = BrooksCorey{FT}();
            _sh2  = VanGenuchten{FT}();

            # Test the struct
            for data_set in [ leaf, root, stem, grass, palm, tree, treet, _vc1, _vc2, _sh1, _sh2]
                @test PkgUtility.FT_test(data_set, FT)
                @test PkgUtility.NaN_test(data_set)
            end;

            # Test Soil types
            _nams = ["Sand", "Loamy Sand", "Sandy Loam", "Loam", "Sandy Clay Loam",
                     "Silt Loam", "Silt", "Clay Loam", "Silty Clay Loam",
                     "Sandy Clay", "Silty Clay", "Clay"];
            for _name in _nams
                create_soil_VC(_sh1, _name);
                create_soil_VC(_sh2, _name);
            end;
        end;
    end;

    @testset "Soil VC" begin
        for FT in [Float32, Float64]
            _sh1  = BrooksCorey{FT}();
            _sh2  = VanGenuchten{FT}();

            # test soil functions
            _rwc = FT(0.5);
            _p   = FT(-0.5);
            for result in [ soil_erwc(_sh2, FT(1)),
                            soil_rwc(_sh1, _p),
                            soil_rwc(_sh2, _p),
                            soil_k_ratio_rwc(_sh1, _rwc),
                            soil_k_ratio_rwc(_sh2, _rwc),
                            soil_k_ratio_swc(_sh1, _rwc*_sh1.Θs),
                            soil_k_ratio_swc(_sh2, _rwc*_sh2.Θs),
                            soil_k_ratio_p25(_sh1, _p),
                            soil_k_ratio_p25(_sh2, _p),
                            soil_p_25_erwc(_sh1, FT(1)),
                            soil_p_25_erwc(_sh2, FT(1)),
                            soil_p_25_rwc(_sh1, _rwc),
                            soil_p_25_rwc(_sh2, _rwc),
                            soil_p_25_swc(_sh1, _rwc*_sh1.Θs),
                            soil_p_25_swc(_sh2, _rwc*_sh2.Θs) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;

            # test if the soil_p_25 and soil_rwc converge
            _p_giv = FT(-2);
            _r_m_1 = soil_erwc(_sh1, _p_giv);
            _r_m_2 = soil_erwc(_sh2, _p_giv);
            _p_m_1 = soil_p_25_erwc(_sh1, _r_m_1);
            _p_m_2 = soil_p_25_erwc(_sh2, _r_m_2);
            @test _p_giv ≈ _p_m_1 ≈ _p_m_2;

            _r_m_1 = soil_rwc(_sh1, _p_giv);
            _r_m_2 = soil_rwc(_sh2, _p_giv);
            _p_m_1 = soil_p_25_rwc(_sh1, _r_m_1);
            _p_m_2 = soil_p_25_rwc(_sh2, _r_m_2);
            @test _p_giv ≈ _p_m_1 ≈ _p_m_2;

            _r_m_1 = soil_swc(_sh1, _p_giv);
            _r_m_2 = soil_swc(_sh2, _p_giv);
            _p_m_1 = soil_p_25_swc(_sh1, _r_m_1);
            _p_m_2 = soil_p_25_swc(_sh2, _r_m_2);
            @test _p_giv ≈ _p_m_1 ≈ _p_m_2;
        end;
    end;

    @testset "Xylem VC" begin
        for FT in [Float32, Float64]
            _vc1 = WeibullSingle{FT}();
            _vc2 = WeibullDual{FT}();
            _sh1 = BrooksCorey{FT}();
            _sh2 = VanGenuchten{FT}();
            _p   = FT(-0.5);
            _v   = FT(0.9);

            # Test xylem_p_crit
            f_st = FT(1)
            for result in [ xylem_p_crit(_vc1, f_st),
                            xylem_p_crit(_vc2, f_st) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;

            # test the xylem_k_ratio
            for result in [ xylem_k_ratio(_vc1, _p, _v),
                            xylem_k_ratio(_vc2, _p, _v),
                            xylem_k_ratio(_vc1, FT(1), _v),
                            xylem_k_ratio(_vc2, FT(1), _v) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Leaf" begin
        for FT in [Float32, Float64]
            leaf  = LeafHydraulics{FT}();
            _f    = FT(0.01);

            # Test the struct
            _lr = xylem_risk(leaf, _f);
            _ec = critical_flow(leaf);

            for result in [_lr, _ec]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Legacy" begin
        for FT in [Float32, Float64]
            _ps   = FT(-0.1);
            _ft   = FT(1);
            _fsl  = FT(0.8);
            _fsh  = FT(0.2);
            _rsl  = FT(0.5);

            # Test the struct
            grass = create_grass(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
            palm  = create_palm(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            tree  = create_tree(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            treet = TreeSimple{FT}();
            pressure_profile!(treet, _ps, _ft);
            @test PkgUtility.FT_test(treet, FT);
            @test PkgUtility.NaN_test(treet);

            treet = TreeSimple{FT}();
            pressure_profile!(treet, _ps, _fsl, _fsh, _rsl);
            @test PkgUtility.FT_test(treet, FT);
            @test PkgUtility.NaN_test(treet);

            for _plant in [grass, palm, tree]
                flow_profile!(_plant);
                pressure_profile!(_plant, SteadyStateMode(); update=true);
                @test PkgUtility.FT_test(_plant, FT);
                @test PkgUtility.NaN_test(_plant);
            end;

            inititialize_legacy!(grass);
            inititialize_legacy!(palm);
            inititialize_legacy!(tree);
            inititialize_legacy!(treet);
            @test PkgUtility.FT_test(treet, FT);
            @test PkgUtility.NaN_test(treet);
        end;
    end;

    @testset "Temperature functions" begin
        for FT in [Float32, Float64]
            leaf  = LeafHydraulics{FT}();
            root  = RootHydraulics{FT}();
            treet = TreeSimple{FT}();
            grass = create_grass(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
            palm  = create_palm(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            tree  = create_tree(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            T     = rand(FT) + 298;

            # test the temperature functions
            temperature_effects!(treet);
            temperature_effects!(grass);
            temperature_effects!(palm );
            temperature_effects!(tree );
            temperature_effects!(leaf, T);
            temperature_effects!(root, T);

            for dataset in [leaf, root, treet, grass, palm, tree]
                @test PkgUtility.FT_test(dataset, FT);
                @test PkgUtility.NaN_test(dataset);
            end;
        end;
    end;

    @testset "Root-related functions" begin
        for FT in [Float32, Float64]
            root  = RootHydraulics{FT}();
            grass = create_grass(FT(-2.1), FT(0.5), collect(FT,0:-0.5:-3.0), collect(FT,0:1:20));

            # test the root q functions
            _p1 = FT(0)
            _p2 = FT(-0.5)
            _p3 = FT(-1.0)
            for result in [ xylem_flow(root, _p1),
                            xylem_flow(root, _p2),
                            xylem_flow(root, _p3) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;

            # test the root qs from q functions
            _f0 = FT(0);
            _f1 = FT(0.1);
            _f2 = FT(0.5);
            _ks = zeros(FT, 5);
            _ps = zeros(FT, 5);
            _qs = zeros(FT, 5);
            for _f in [_f0, _f1, _f2]
                roots_flow!(grass.roots, _ks, _ps, _qs, _f);
                @test PkgUtility.FT_test(grass, FT);
                @test PkgUtility.NaN_test(grass);
            end;

            # test the recalculate_roots_flow function
            _ks = zeros(FT, 5);
            _ps = zeros(FT, 5);
            _qs = zeros(FT, 5);
            roots_flow!(grass.roots, _ks, _ps, _qs, FT(0.5));
            roots_flow!(grass, FT(0.5));
            @test PkgUtility.NaN_test(grass);
        end;
    end;

    @testset "Pressure functions" begin
        for FT in [Float32, Float64]
            leaf  = LeafHydraulics{FT}();
            root  = RootHydraulics{FT}();
            stem  = StemHydraulics{FT}();
            grass = create_grass(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
            palm  = create_palm(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            tree  = create_tree(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            treet = TreeSimple{FT}();

            # test the end_pressure function
            _f_1 = FT(0.001)
            _f_2 = FT(1)
            for result in [ end_pressure(leaf, _f_1),
                            end_pressure(leaf, _f_2),
                            end_pressure(root, _f_1),
                            end_pressure(root, _f_2),
                            end_pressure(stem, _f_1),
                            end_pressure(stem, _f_2) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;

            # test the end_pressure function for plants
            for result in [ end_pressure(treet, _f_1),
                            end_pressure(treet, _f_2),
                            end_pressure(treet, _f_1, _f_1, FT(0.5)) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;
        end;
    end;

    @testset "Plant level" begin
        for FT in [Float32, Float64]
            grass = create_grass(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
            palm  = create_palm(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            tree  = create_tree(FT(-2.1), FT(5.5), FT(6), FT[0,-1,-2,-3], collect(FT,0:1:20));
            treet = TreeSimple{FT}();

            # test the critical_flow function
            for result in [ critical_flow(treet) ]
                @test PkgUtility.FT_test(result, FT);
                @test PkgUtility.NaN_test(result);
            end;

            # test the plant conductances
            plant_conductances!(treet);
        end;
    end;

    @testset "Capacitance" begin
        for FT in [Float32, Float64]
            grass = create_grass(FT(-2.1), FT(0.5), FT[0,-1,-2,-3], collect(FT,0:1:20));
            palm  = create_palm(FT(-2.1), FT(5.4), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));
            tree  = create_tree(FT(-2.1), FT(5.4), FT(8), FT[0,-1,-2,-3], collect(FT,0:1:20));

            # test the critical_flow function
            update_PVF!(grass, FT(1)); @test true;
            update_PVF!(palm , FT(1)); @test true;
            update_PVF!(tree , FT(1)); @test true;
            for leaf in tree.leaves
                leaf.q_out = 2e-3;
            end
            update_PVF!(tree , FT(1)); @test true;
            update_PVF!(tree.roots[1], FT(1e6));
            update_PVF!(tree.roots[1], FT(1e6), true);
            update_PVF!(tree.trunk, FT(1e6));
            update_PVF!(tree.leaves[1], FT(1e6));

            # extra tests
            p_from_volume(PVCurveSegmented{FT}(), FT(0.9), FT(300));
            p_from_volume(PVCurveSegmented{FT}(), FT(0.6), FT(300));
            p_from_volume(PVCurveSegmented{FT}(), FT(0.1), FT(300));
        end;
    end;
end;
