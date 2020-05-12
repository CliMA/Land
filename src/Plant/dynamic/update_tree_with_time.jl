# function to update the tree with time, default at a Δt=1s using Wang-2020 scheme
function update_tree_with_time(tree::StructTree, Δt::Number=1.0; scheme::String="Wang2020")
    # 1. use the gsw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi::StructTreeCanopyLayer in tree.canopy.canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        for leaf in canopyi.leaf_list
            d_leaf       = (get_saturated_vapor_pressure(leaf.t_leaf) - canopyi.p_h2o) / canopyi.p_atm
            leaf.e       = leaf.gsw * d_leaf
            anagrpi      = get_leaf_an_ag_r_pi_from_gsc(
                                                        v25 = leaf.v_max,
                                                        j25 = leaf.j_max,
                                                     Γ_star = leaf.Γ_star,
                                                        gsc = leaf.gsc,
                                                        p_a = canopyi.p_a,
                                                        tem = leaf.t_leaf,
                                                        par = leaf.par,
                                                      p_atm = canopyi.p_atm,
                                                       p_O₂ = canopyi.p_O₂,
                                                        r25 = leaf.r_25,
                                                       unit = "K")
            leaf.a_net   = anagrpi[1]
            leaf.a_gross = anagrpi[2]
            leaf.r       = anagrpi[3]
            leaf.p_i     = anagrpi[4]
        end

        # 1.2 update the e_list and q_list for the canopy layer, per tree
        canopyi.e_list = [leaf.e for leaf in canopyi.leaf_list]
        canopyi.q_list = canopyi.e_list .* canopyi.la_list
    end

    # 2 update the pressure profiles for roots
    q_canopy_list      = [sum(canopyi.q_list) for canopyi in tree.canopy.canopy_list]
    q_sum              = sum(q_canopy_list)
    p_base,q_root_list = get_p_base_q_list_from_q(tree, q_sum)
    for indx in 1:length(tree.roots.root_list)
        rooti = tree.roots.root_list[indx]
        update_struct_from_q(rooti, q_root_list[indx])
    end

    # 3 update the pressure profile in the trunk
    tree.trunk.p_ups = p_base
    update_struct_from_q(tree.trunk, q_sum)

    # 4 update the pressure profiles in the branches and leaves
    for indx in 1:length(tree.branch.branch_list)
        # 4.1 update the pressure profile in the branch
        branchi       = tree.branch.branch_list[indx]
        canopyi       = tree.canopy.canopy_list[indx]
        branchi.p_ups = tree.trunk.p_dos
        update_struct_from_q(branchi, q_canopy_list[indx])
        
        # 4.2 update the pressure profiles for leaves
        for leaf in canopyi.leaf_list
            leaf.p_ups = branchi.p_dos
            update_struct_from_q(leaf, leaf.e)
        end
    end

    # 5. determine how much gsw and gsc should change with time, use Wang 2020 model for day time, will add scheme for nighttime later
    for canopyi in tree.canopy.canopy_list
        for indx in 1:length(canopyi.leaf_list)
            leaf = canopyi[indx]
            ∂A∂E = get_marginal_gain(canopyi, leaf)
            if scheme=="Wang2020"
                ∂Θ∂E = get_marginal_penalty_wang(leaf)
            else # always maximize A otherwise
                ∂Θ∂E = 0.0
            end
            leaf.gsw += leaf.gs_nssf * (∂A∂E - ∂Θ∂E) * Δt
            leaf.gsc  = leaf.gsw / 1.6 / ( 1.0 + leaf.g_ias_c * leaf.gsw^(leaf.g_ias_e) )
        end
    end
end