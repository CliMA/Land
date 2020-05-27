"""
    update_tree_with_time!(tree::Tree, Δt::FT, scheme::AbstractOptimizationStomatalModel; updating::Bool=false)

Update the tree flux information and hydraulic parameters, given
- `tree` One [`Tree`](@ref) type
- `Δt` Δt of non-steady state gs
- `scheme` An `AbstractOptimizationStomatalModel` type optimization stomatal model scheme
- `updating` If true, plant hydraulic parameters will be updated for the [`Tree`](@ref)

The function updates the stomatal conductance via a `Δgsw = factor * (∂A/∂E - ∂Θ/∂E)` for optimization stomatal models. The `∂A/∂E` is the same for every stomatal control model, but `∂Θ/∂E` differs among models.

"""
function update_tree_with_time!(tree::Tree, Δt::FT, scheme::AbstractOptimizationStomatalModel; updating::Bool=false) where {FT}
    # 0. unpack necessary structs
    @unpack branch_list,canopy_list,root_list   = tree

    # 1. use the gsw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi in canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        canopyi.d_list  = (get_saturated_vapor_pressure(canopyi.t_list) .- canopyi.p_H₂O) ./ canopyi.p_atm
        canopyi.e_list  = canopyi.gsw_list .* canopyi.d_list
        canopyi.q_list  = canopyi.e_list   .* canopyi.la_list
        anagrpi_lists   = get_an_ag_r_pi_from_gsc_list(
                                                       tree.photo_para_set,
                                                       canopyi.gsc_list,
                                                       canopyi.v_max,
                                                       canopyi.j_max,
                                                       canopyi.p_max,
                                                       canopyi.p_a,
                                                       canopyi.t_list,
                                                       canopyi.par_list,
                                                       canopyi.p_atm,
                                                       canopyi.p_O₂,
                                                       canopyi.r_25,
                                                       canopyi.curvature,
                                                       canopyi.qy)
        canopyi.an_list = anagrpi_lists[1]
        canopyi.ag_list = anagrpi_lists[2]
        canopyi.r_list  = anagrpi_lists[3]
        canopyi.pi_list = anagrpi_lists[4]
    end

    # update the pressure profile in the trunk, branch, and leaf
    if updating
        # 2. update the pressure profiles for roots
        q_canopy_list      = [sum(canopyi.q_list) for canopyi in canopy_list]
        q_sum              = sum(q_canopy_list)
        p_base,q_root_list = get_p_base_q_list_from_q(tree, q_sum)
        for indx in 1:length(root_list)
            rooti = root_list[indx]
            update_struct_from_q!(rooti, q_root_list[indx])
        end

        # 3. update the pressure profile in the trunk
        tree.trunk.p_ups = p_base
        update_struct_from_q!(tree.trunk, q_sum)

        # 4. update the pressure profiles in the branches and leaves
        for indx in 1:length(branch_list)
            # 4.1 update the pressure profile in the branch
            branchi       = branch_list[indx]
            canopyi       = canopy_list[indx]
            branchi.p_ups = tree.trunk.p_dos
            update_struct_from_q!(branchi, q_canopy_list[indx])
            
            # 4.2 update the pressure profiles for leaves
            for ind_leaf in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[ind_leaf]
                leaf.p_ups = branchi.p_dos
                update_struct_from_q!(leaf, canopyi.e_list[ind_leaf])
            end
        end
    end

    # 5. determine how much gsw and gsc should change with time
    for canopyi in canopy_list
        # 5.1 compute the ∂A∂E and ∂Θ∂E for each leaf
        list_∂A∂E = get_marginal_gain(canopyi, tree.photo_para_set)
        list_∂Θ∂E = get_marginal_penalty(canopyi, tree.stomata_scheme)

        # update the gsw and gsc for each leaf
        canopyi.gsw_list += canopyi.gs_nssf .* (list_∂A∂E - list_∂Θ∂E) .* Δt
        canopyi.gsw_list[ canopyi.gsw_list .< 0 ] .= 0
        canopyi.gsc_list  = canopyi.gsw_list ./ FT(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.gsw_list .^ (canopyi.g_ias_e) )
    end
end




"""
    update_tree_with_time!(tree::Tree, Δt::FT, scheme::AbstractEmpiricalStomatalModel; updating::Bool=false)

Update the tree flux information and hydraulic parameters, given
- `tree` One [`Tree`](@ref) type
- `Δt` Δt of non-steady state gs
- `scheme` An `AbstractEmpiricalStomatalModel` type empirical stomatal model scheme
- `updating` If true, plant hydraulic parameters will be updated for the [`Tree`](@ref)

The function updates the stomatal conductance via a `Δgsw = factor * (gs_ss - gs_nss)` for empirical stomatal models like Ball-Berry, Leuning, Medlyn, and Gentine models. The `gs_ss` stands for the gs at steady state, and `gs_nss` stands for gs at non-steady state.

Need to link to empirical models before use it.

"""
function update_tree_with_time!(tree::Tree, Δt::FT, scheme::AbstractEmpiricalStomatalModel; updating::Bool=false) where {FT}
    # 0. unpack necessary structs
    @unpack branch_list,canopy_list,root_list   = tree

    # 1. use the gsw from last time instant and update e and a first, because Tleaf was from the last instant
    for canopyi in canopy_list
        # 1.1 update the e and a for each leaf, per leaf area
        canopyi.d_list  = (get_saturated_vapor_pressure(canopyi.t_list) .- canopyi.p_H₂O) ./ canopyi.p_atm
        canopyi.e_list  = canopyi.gsw_list .* canopyi.d_list
        canopyi.q_list  = canopyi.e_list   .* canopyi.la_list
        anagrpi_lists   = get_an_ag_r_pi_from_gsc_list(
                                                       tree.photo_para_set,
                                                       canopyi.gsc_list,
                                                       canopyi.v_max,
                                                       canopyi.j_max,
                                                       canopyi.p_max,
                                                       canopyi.p_a,
                                                       canopyi.t_list,
                                                       canopyi.par_list,
                                                       canopyi.p_atm,
                                                       canopyi.p_O₂,
                                                       canopyi.r_25,
                                                       canopyi.curvature,
                                                       canopyi.qy)
        canopyi.an_list = anagrpi_lists[1]
        canopyi.ag_list = anagrpi_lists[2]
        canopyi.r_list  = anagrpi_lists[3]
        canopyi.pi_list = anagrpi_lists[4]
    end

    # update the pressure profile in the trunk, branch, and leaf
    if updating
        # 2. update the pressure profiles for roots
        q_canopy_list      = [sum(canopyi.q_list) for canopyi in canopy_list]
        q_sum              = sum(q_canopy_list)
        p_base,q_root_list = get_p_base_q_list_from_q(tree, q_sum)
        for indx in 1:length(root_list)
            rooti = root_list[indx]
            update_struct_from_q!(rooti, q_root_list[indx])
        end

        # 3. update the pressure profile in the trunk
        tree.trunk.p_ups = p_base
        update_struct_from_q!(tree.trunk, q_sum)

        # 4. update the pressure profiles in the branches and leaves
        for indx in 1:length(branch_list)
            # 4.1 update the pressure profile in the branch
            branchi       = branch_list[indx]
            canopyi       = canopy_list[indx]
            branchi.p_ups = tree.trunk.p_dos
            update_struct_from_q!(branchi, q_canopy_list[indx])
            
            # 4.2 update the pressure profiles for leaves
            for ind_leaf in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[ind_leaf]
                leaf.p_ups = branchi.p_dos
                update_struct_from_q!(leaf, canopyi.e_list[ind_leaf])
            end
        end
    end

    # 5. determine how much gsw and gsc should change with time
    for canopyi in canopy_list
        # 5.1 compute the gs at steady state, function to be added
        # the gs_nssf value may not work for the empirical models (too big!)
        gss_list = canopyi.gsw_list .* 1

        # update the gsw and gsc for each leaf
        canopyi.gsw_list += canopyi.gs_nssf .* (gss_list - canopyi.gsw_list) .* Δt
        canopyi.gsw_list[ canopyi.gsw_list .< 0 ] .= 0
        canopyi.gsc_list  = canopyi.gsw_list ./ FT(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.gsw_list .^ (canopyi.g_ias_e) )
    end
end
