# function to calculate marginal gain for any optimization model
function get_marginal_gain(canopy_layer::struct_tree_canopy_layer, indx)
    leaf        = canopy_layer[indx]
    leaf_e1     = leaf.e + 1E-6
    leaf_d      = (get_saturated_vapor_pressure(leaf.t_leaf) - canopy_layer.p_h2o) /  canopy_layer.p_atm
    leaf_gsw1   = leaf.e1 / leaf_d
    leaf_gsc1   = leaf_gsw1 / 1.6 / ( 1.0 + leaf.g_ias_c * leaf_gsw1^(leaf.g_ias_e) )
    leaf_a1,p_i = get_leaf_an_pi_from_gsc(leaf.v_max, leaf.j_max, leaf.g_star, leaf_gsc1, canopy_layer.p_a, leaf.t_leaf, leaf.par, canopy_layer.p_atm, canopy_layer.p_o2, leaf.respir)
    dade        = leaf_a1 - leaf.a_net # mol CO2 mol^-1 H2O, 1E-6 mol (H2O) and umol = 1E-6 mol (CO2) cancel out
    return dade
end