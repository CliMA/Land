# function to get marginal carbon penalty for the Wang 2020 model
function get_marginal_penalty_wang(canopy_layer::struct_tree_canopy_layer, indx)
    leaf = canopy_layer[indx]
    dTde = leaf.a_net / (leaf.e_crit - leaf.e) * 1E-6
    return dTde
end