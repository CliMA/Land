# test the structures
@info "Testing FT and NaN for the structures...";
@testset "StomtaModels --- structures" begin
    for FT in [Float32, Float64]
        can_3 = CanopyLayer{FT}(n_leaf=2);
        esm_1 = ESMBallBerry{FT}();
        esm_2 = ESMGentine{FT}();
        esm_3 = ESMLeuning{FT}();
        esm_4 = ESMMedlyn{FT}();
        osm_1 = OSMEller{FT}();
        osm_2 = OSMSperry{FT}();
        osm_3 = OSMWang{FT}();
        osm_4 = OSMWAP{FT}();
        osm_5 = OSMWAPMod{FT}();

        # test the structures
        for result in [ can_3, esm_1, esm_2, esm_3, esm_4, osm_4, osm_5 ]
            @test FT_test(result, FT);
        end
    end
end



# test the structures
println();
@info "Testing gas exchange functions...";
@testset "StomtaModels --- gas exchange functions" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}("C3");
        leaf_4 = Leaf{FT}("C4");
        rand_T = rand(FT) + 298;
        can_3  = CanopyLayer{FT}(n_leaf=2);
        can_4  = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the refresh functions
        update_leaf_TP!(can_3, hs, envir);
        update_leaf_TP!(can_4, hs, envir);
        @test true;

        update_leaf_AK!(can_3, hs, envir);
        update_leaf_AK!(can_4, hs, envir);
        @test true;

        gas_exchange!(can_3, envir, GlcDrive(), 1, FT(0.1));
        gas_exchange!(can_4, envir, GlcDrive(), 1, FT(0.1));
        @test true;

        gas_exchange!(can_3, envir, GswDrive(), 1, FT(0.05));
        gas_exchange!(can_4, envir, GswDrive(), 1, FT(0.05));
        @test true;

        can_3.g_sw[2] = 0;
        can_4.g_sw[2] = 0;
        gsw_control!(can_3, envir, 2);
        gsw_control!(can_4, envir, 2);
        @test true;
    end
end



# test the structures
println();
@info "Testing empirical formulations...";
@testset "StomtaModels --- empirical formulations" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}("C3");
        leaf_4 = Leaf{FT}("C4");
        rand_T = rand(FT) + 298;
        can_3  = CanopyLayer{FT}(n_leaf=2);
        can_4  = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();

        # test the empirical model formulations
        for result in [ stomatal_conductance(esm_1, leaf_3, envir, FT(1)),
                        stomatal_conductance(esm_2, leaf_3, envir, FT(1)),
                        stomatal_conductance(esm_3, leaf_3, envir, FT(1)),
                        stomatal_conductance(esm_4, leaf_3, envir, FT(1)),
                        stomatal_conductance(esm_1, can_3, envir, FT(1)),
                        stomatal_conductance(esm_2, can_3, envir, FT(1)),
                        stomatal_conductance(esm_3, can_3, envir, FT(1)),
                        stomatal_conductance(esm_4, can_3, envir, FT(1)),
                        stomatal_conductance(esm_1, can_3, envir, FT(1), 1),
                        stomatal_conductance(esm_2, can_3, envir, FT(1), 1),
                        stomatal_conductance(esm_3, can_3, envir, FT(1), 1),
                        stomatal_conductance(esm_4, can_3, envir, FT(1), 1) ]
            @test FT_test(result, FT);
            @test NaN_test(result);
        end
    end
end



# test the structures
println();
@info "Testing stomatal models...";
@testset "StomtaModels --- stomatal models" begin
    for FT in [Float32, Float64]
        envir  = AirLayer{FT}();
        leaf_3 = Leaf{FT}("C3");
        leaf_4 = Leaf{FT}("C4");
        rand_T = rand(FT) + 298;
        can_3  = CanopyLayer{FT}(n_leaf=2);
        can_4  = CanopyLayer{FT}(n_leaf=2);
        esm_1  = ESMBallBerry{FT}();
        esm_2  = ESMGentine{FT}();
        esm_3  = ESMLeuning{FT}();
        esm_4  = ESMMedlyn{FT}();
        osm_1  = OSMEller{FT}();
        osm_2  = OSMSperry{FT}();
        osm_3  = OSMWang{FT}();
        osm_4  = OSMWAP{FT}();
        osm_5  = OSMWAPMod{FT}();
        hs     = LeafHydraulics{FT}();
        svc    = VanGenuchten{FT}();

        # test the solution functions
        for can in [can_3, can_4]
            for sm in [esm_1, esm_2, esm_3, esm_4]
                for result in [ solution_diff!(FT(0.1), can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}(), GlcDrive(), 1),
                                solution_diff!(FT(0.1), can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}(), GlcDrive(), 1),
                                solution_diff!(FT(0.1), can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}(), GlcDrive(), 1),
                                solution_diff!(FT(0.1), can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}(), GlcDrive(), 1),
                                solution_diff!(FT(0.1), can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}(), GlcDrive(), 1),
                                solution_diff!(FT(0.1), can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}(), GlcDrive(), 1),
                                ]
                    @test FT_test(result, FT);
                    @test NaN_test(result);
                end
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                result = solution_diff!(FT(0.1), can, hs, envir, sm, GlcDrive(), 1);
                @test FT_test(result, FT);
                @test NaN_test(result);
            end
        end

        # test the stomata solutions
        for can in [can_3, can_4]
            for sm in [esm_1, esm_2, esm_3, esm_4]
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearKleaf{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearKsoil{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearKleaf{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearKsoil{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}(), 1);
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}(), 1);
                @test true;
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                gas_exchange!(can, hs, envir, sm, 1);
                @test true;
            end
            gas_exchange!(can, TreeSimple{FT}(), envir, osm_3);
            @test true;
        end

        for can in [can_3, can_4]
            for sm in [esm_1, esm_2, esm_3, esm_4]
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearKleaf{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearKsoil{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearPleaf{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearPsoil{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaGLinearSWC{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearKleaf{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearKsoil{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearPleaf{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearPsoil{FT}());
                @test true;
                gas_exchange!(can, hs, svc, FT(-1), FT(0.4), envir, sm, BetaVLinearSWC{FT}());
                @test true;
            end
            for sm in [osm_1, osm_2, osm_3, osm_4, osm_5]
                gas_exchange!(can, hs, envir, sm);
                @test true;
            end
        end

        # test the nighttime stomatal conductance solution for Wang model
        can_3.APAR[1] = 0;
        gas_exchange!(can_3, hs, envir, osm_3);

        # test the prognostic g_sw functions
        prognostic_gsw!(can_3, envir, esm_1, FT(1), FT(120));
        prognostic_gsw!(can_3, hs, envir, osm_3, FT(120));
        @test true;
    end
end
