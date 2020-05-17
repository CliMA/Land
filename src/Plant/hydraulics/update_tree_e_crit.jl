using DocStringExtensions

"""
    update_tree_e_crit!(tree, displaying)
This function update e_crit for each layer (sunlit and shaded).
The e_crit for each leaf is either increased or decreased so that p_leaf is in the range -Inf to -20.0 MPa.
The e_crit increases by 1E-6 mol m^-2 s^-1 if p_leaf is less negative than -20 MPa.
The e_crit decreases by 1E-7 mol m^-2 s^-1 if p_leaf is -Inf.
There is a "displaying" variable to control whether the information of q_sum will display (true) or not (false).
"""
function update_tree_e_crit!(tree::Tree; displaying::Bool=false)
    # unpack necessary structs
    @unpack branch_list = tree.branch
    @unpack canopy_list = tree.canopy
    @unpack root_list   = tree.roots

    while true
        # define judge = 0
        judge = 0

        # q for the whole tree is the sum of q in each leaves
        q_branch_list = [sum(canopyi.ec_list .* canopyi.la_list) for canopyi in canopy_list]
        q_sum = sum( q_branch_list )

        # calculate the p_base from q_sum
        p_base,q_list = get_p_base_q_list_from_q(tree, q_sum)

        # calculate the trunk-branch joint pressure from q_sum and p_base
        p_trunk_branch = get_struct_p_end_from_q(tree.trunk, q_sum; p_ini=p_base)

        # for each canopy layer
        for indx in 1:length(branch_list)
            branchi = branch_list[indx]
            canopyi = canopy_list[indx]
            # compute branch leaf joint pressure for each canopy layer
            p_branch_leaf = get_struct_p_end_from_q(branchi, q_branch_list[indx]; p_ini=p_trunk_branch)

            # compute p_leaf for each leaf
            for indy in 1:length(canopyi.leaf_list)
                leaf = canopyi.leaf_list[indy]
                p_leaf = get_struct_p_end_from_q(leaf, canopyi.ec_list[indy]; p_ini=p_branch_leaf)
                # increas e_crit by 1E-6 or decrease it by 1E-7
                if p_leaf >= -20.0
                    judge += 1
                    canopyi.ec_list[indy] += ΔEP_6
                elseif p_leaf == -Inf
                    judge += 1
                    canopyi.ec_list[indy] -= ΔEP_7
                end
            end
        end

        # if judge > 0, continue
        if judge==0
            break
        else
            if displaying
                println("The Q_sum was ", q_sum, " and then ", judge, " modification(s) was(were) made.")
            end
        end
    end
end