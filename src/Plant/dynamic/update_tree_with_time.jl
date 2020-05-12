function update_tree_with_time(tree, dtime=1.0)
    # 1. use the gsw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi in tree.canopy.canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        for leaf in canopyi.leaf_list
            d_leaf       = (get_saturated_vapor_pressure(leaf.t_leaf) - canopyi.p_h2o) / canopyi.p_atm
            leaf.e       = leaf.gsw * d_leaf
            anpiag       = get_leaf_an_pi_ag_from_gsc(leaf.v_max, leaf.j_max, leaf.g_star, leaf.gsc, canopyi.p_a, leaf.t_leaf, leaf.par, canopyi.p_atm, canopyi.p_o2, leaf.respir)
            leaf.a_net   = anpiag[1]
            leaf.p_i     = anpiag[2]
            leaf.a_gross = anpiag[3]
            leaf.r       = anpiag[3] - anpiag[1]
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
        update_struct_from_q_rhizosphere(rooti, q_root_list[indx])
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
end