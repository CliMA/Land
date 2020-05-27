using Land.Plant

PLT = Plant

@testset "Plant - isNaN and FT Consistency" begin
    for FT in [Float32 Float64]
        tree = PLT.Tree{FT,5,20,325}();
           v25::FT = FT(80.0)
           j25::FT = FT(135.0)
           p25::FT = FT(120.0)
           gsc::FT = FT(0.1)
           p_a::FT = FT(40.0)
           p_i::FT = FT(30.0)
           par::FT = FT(1000.0)
           tem::FT = FT(298.15)
         p_atm::FT = FT(101325.0)
          p_O₂::FT = FT(21278.25)
           r25::FT = FT(Inf)
     curvature::FT = FT(0.9)
            qy::FT= FT(0.42)

        # test the FT Consistency
        for tmp in [PLT.get_p_base_q_list_from_q(tree, FT(1.0)),
                    PLT.get_q_layer_from_p_base(tree.root_list[1], FT(-0.1)),
                    PLT.get_struct_p_end_from_q(tree.root_list[1], FT(0.1); p_ini=FT(Inf)),
                    PLT.get_struct_p_end_from_q(tree.trunk, FT(1.0); p_ini=FT(Inf)),
                    PLT.get_struct_p_end_from_q(tree.branch_list[1], FT(0.1); p_ini=FT(Inf)),
                    PLT.get_struct_p_end_from_q(tree.canopy_list[1].leaf_list[1], FT(0.05); p_ini=FT(Inf)),
                    PLT.get_marginal_gain(tree.canopy_list[1], 1, tree.photo_para_set),
                    PLT.get_marginal_gain(tree.canopy_list[1], tree.photo_para_set),
                    PLT.get_marginal_penalty(tree.canopy_list[1], 1, tree.stomata_scheme),
                    PLT.get_marginal_penalty(tree.canopy_list[1], tree.stomata_scheme),
                    PLT.get_a_par_curve(tree.photo_para_set, v25, j25, p25, gsc, p_a, tem, p_atm, p_O₂, r25, curvature, qy),
                    PLT.get_a_pi_curve(tree.photo_para_set, v25, j25, p25, tem, par, p_O₂,r25, curvature, qy),
                    PLT.get_relative_surface_tension(FT(298.0)),
                    PLT.get_relative_viscosity(FT(298.0)),
                    PLT.get_saturated_vapor_pressure(FT(298.0)),
                    PLT.get_saturated_vapor_pressure(ones(FT,10) .* FT(298.0)),
                    PLT.get_specific_latent_heat(FT(298.5))]
            if typeof(tmp) <: Number
                @test typeof(tmp) == FT
                @test !isnan(tmp)
            else
                for i in tmp
                    if typeof(i) <: Number
                        @test typeof(i) == FT
                        @test !isnan(i)
                    else
                        for j in i
                            if typeof(j) <: Number
                                @test typeof(j) == FT
                                @test !isnan(j)
                            end
                        end
                    end
                end
            end
        end
    end
end