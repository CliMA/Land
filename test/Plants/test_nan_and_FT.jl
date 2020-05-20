using Land.Plant

@testset "Plant - isNaN and FT Consistency" begin
    for FT in [Float32 Float64]
        tree = Plant.Tree{FT,5,20,325}();
           v25::FT = FT(80.0)
           j25::FT = FT(135.0)
        Γ_star::FT = FT(2.5)
           gsc::FT = FT(0.1)
           p_a::FT = FT(40.0)
           p_i::FT = FT(30.0)
           par::FT = FT(1000.0)
           tem::FT = FT(298.15)
         p_atm::FT = FT(101325.0)
          p_O₂::FT = FT(21278.25)
           r25::FT = FT(Inf)

        # test the FT Consistency
        p_base,q_list = Plant.get_p_base_q_list_from_q(tree, FT(1.0))
        @test typeof(p_base) == FT;
        @test !isnan(p_base);
        for tmp in q_list
            @test typeof(tmp) == FT;
            @test !isnan(tmp);
        end

        q_layer = Plant.get_q_layer_from_p_base(tree.roots.root_list[1], FT(-0.1));
        @test typeof(q_layer) == FT;
        @test !isnan(q_layer);

        p_end = Plant.get_struct_p_end_from_q(tree.roots.root_list[1], FT(0.1); p_ini=FT(Inf));
        @test typeof(p_end) == FT;
        @test !isnan(p_end);

        p_end = Plant.get_struct_p_end_from_q(tree.trunk, FT(1.0); p_ini=FT(Inf));
        @test typeof(p_end) == FT;
        @test !isnan(p_end);

        p_end = Plant.get_struct_p_end_from_q(tree.branch.branch_list[1], FT(0.1); p_ini=FT(Inf));
        @test typeof(p_end) == FT;
        @test !isnan(p_end);

        p_end = Plant.get_struct_p_end_from_q(tree.canopy.canopy_list[1].leaf_list[1], FT(0.05); p_ini=FT(Inf));
        @test typeof(p_end) == FT;
        @test !isnan(p_end);

        temp = Plant.get_leaf_an_ag_r_from_pi(
                                              v25 = v25,
                                              j25 = j25,
                                           Γ_star = Γ_star,
                                              p_i = p_i,
                                              tem = tem,
                                              par = par,
                                             p_O₂ = p_O₂,
                                              r25 = r25);
        for tmp in temp
            @test typeof(tmp) == FT;
            @test !isnan(tmp);
        end

        temp = Plant.get_leaf_an_ag_r_pi_from_gsc(
                                                  v25 = v25,
                                                  j25 = j25,
                                               Γ_star = Γ_star,
                                                  gsc = gsc,
                                                  p_a = p_a,
                                                  tem = tem,
                                                  par = par,
                                                p_atm = p_atm,
                                                 p_O₂ = p_O₂,
                                                  r25 = r25);
        for tmp in temp
            @test typeof(tmp) == FT;
            @test !isnan(tmp);
        end

        temp = Plant.get_leaf_an_ag_r_pi_from_gsc_list(
                                                       v25 = v25,
                                                       j25 = j25,
                                                    Γ_star = Γ_star,
                                                  gsc_list = FT(0.1) .* ones(FT,10),
                                                       p_a = p_a,
                                                  tem_list = FT(298.15) .* ones(FT,10),
                                                  par_list = FT(1000.0) .* ones(FT,10),
                                                     p_atm = p_atm,
                                                      p_O₂ = p_O₂,
                                                       r25 = r25);
        for tmp1 in temp
            for tmp2 in tmp1
                @test typeof(tmp2) == FT;
                @test !isnan(tmp2);
            end
        end

        temp = Plant.get_leaf_j(FT(100.0), FT(1000.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_leaf_jmax(FT(100.0), FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_leaf_r_from_r25(FT(2.0), FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_leaf_r_from_v25(FT(100.0), FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_leaf_vcmax(FT(100.0), FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_marginal_gain(tree.canopy.canopy_list[1], 1);
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_marginal_gain(tree.canopy.canopy_list[1]);
        for tmp in temp
            @test typeof(tmp) == FT;
            @test !isnan(tmp);
        end

        temp = Plant.get_marginal_penalty_wang(tree.canopy.canopy_list[1], 1);
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_marginal_penalty_wang(tree.canopy.canopy_list[1]);
        for tmp in temp
            @test typeof(tmp) == FT;
            @test !isnan(tmp);
        end

        temp = Plant.get_a_par_curve(
                                     v25 = v25,
                                     j25 = j25,
                                  Γ_star = Γ_star,
                                     gsc = gsc,
                                     p_a = p_a,
                                     tem = tem,
                                   p_atm = p_atm,
                                    p_O₂ = p_O₂,
                                     r25 = r25);
        for tmp1 in temp
            for tmp2 in tmp1
                @test typeof(tmp2) == FT;
                @test !isnan(tmp2);
            end
        end

        temp = Plant.get_a_pi_curve(;
                                    v25 = v25,
                                    j25 = j25,
                                 Γ_star = Γ_star,
                                    tem = tem,
                                    par = par,
                                   p_O₂ = p_O₂,
                                    r25 = r25);
        for tmp1 in temp
            for tmp2 in tmp1
                @test typeof(tmp2) == FT;
                @test !isnan(tmp2);
            end
        end

        temp = Plant.get_relative_surface_tension(FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_relative_viscosity(FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_saturated_vapor_pressure(FT(298.0));
        @test typeof(temp) == FT;
        @test !isnan(temp);

        temp = Plant.get_saturated_vapor_pressure(ones(FT,10) .* FT(298.0));
        for tmp in temp
            @test typeof(tmp) == FT;
            @test !isnan(tmp);
        end

        temp = Plant.get_specific_latent_heat(FT(298.5));
        @test typeof(temp) == FT;
        @test !isnan(temp);
    end
end