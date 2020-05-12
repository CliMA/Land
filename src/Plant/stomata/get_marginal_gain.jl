# function to calculate marginal gain for any optimization model
function get_marginal_gain(canopyi::StructTreeCanopyLayer, leaf::StructTreeLeaf)
    leaf_e1   = leaf.e + 1E-6
    leaf_d    = (get_saturated_vapor_pressure(leaf.t_leaf) - canopyi.p_H₂O) /  canopyi.p_atm
    leaf_gsw1 = leaf_e1 / leaf_d
    leaf_gsc1 = leaf_gsw1 / 1.6 / ( 1.0 + leaf.g_ias_c * leaf_gsw1^(leaf.g_ias_e) )
    anagrpi   = get_leaf_an_ag_r_pi_from_gsc(
                                             v25 = leaf.v_max,
                                             j25 = leaf.j_max,
                                          Γ_star = leaf.Γ_star,
                                             gsc = leaf_gsc1,
                                             p_a = canopyi.p_a,
                                             tem = leaf.t_leaf,
                                             par = leaf.par,
                                           p_atm = canopyi.p_atm,
                                            p_O₂ = canopyi.p_O₂,
                                             r25 = leaf.r_25,
                                            unit = "K")
    leaf_a1   = anagrpi[1]
    ∂A∂E      = leaf_a1 - leaf.a_net # mol CO₂ mol⁻¹ H₂O, 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return ∂A∂E
end