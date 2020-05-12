# function to update e_crit for each layer (sunlit and shaded), increase or decrease by 1E-6 mol m^-2 s^-1
function update_tree_e_crit(tree::StructTree; displaying::Bool=false)
    while true
        # define judge = 0
        judge = 0

        # q for the whole tree is the sum of q in each leaves
        q_branch_list = []
        for indx in 1:length(tree.branch.branch_list)
            canopyi = tree.canopy.canopy_list[indx]
            e_crit_list  = [leaf.e_crit for leaf in canopyi.leaf_list]
            push!( q_branch_list, sum( e_crit_list .* canopyi.la_list) )
        end
        q_sum = sum( q_branch_list )

        # calculate the p_base from q_sum
        p_base,q_list = get_p_base_q_list_from_q(tree, q_sum)

        # calculate the trunk-branch joint pressure from q_sum and p_base
        p_trunk_branch = get_struct_p_end_from_q(tree.trunk, q_sum; p_ini=p_base)

        # for each canopy layer
        for indx in 1:length(tree.branch.branch_list)
            branchi = tree.branch.branch_list[indx]
            canopyi = tree.canopy.canopy_list[indx]
            # compute branch leaf joint pressure for each canopy layer
            p_branch_leaf = get_struct_p_end_from_q(branchi, q_branch_list[indx]; p_ini=p_trunk_branch)

            # compute p_leaf for each leaf
            for leaf in canopyi.leaf_list
                p_leaf = get_struct_p_end_from_q(leaf, leaf.e_crit; p_ini=p_branch_leaf)
                # increas e_crit by 1E-6 or decrease it by 1E-7
                if p_leaf >= -20.0
                    judge += 1
                    leaf.e_crit += 1E-6
                elseif p_leaf == -Inf
                    judge += 1
                    leaf.e_crit -= 1E-7
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