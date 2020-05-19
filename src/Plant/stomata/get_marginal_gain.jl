"""
    get_marginal_gain(canopyi, indx)

Marginal water use efficiency `∂A/∂E` for a leaf in a canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
"""
function get_marginal_gain(canopyi::CanopyLayer, indx::Number)
    FTYP      = eltype(canopyi.an_list[indx])
    leaf_e1   = canopyi.e_list[indx] + FTYP(1E-6)
    leaf_d    = (get_saturated_vapor_pressure(canopyi.t_list[indx]) - canopyi.p_H₂O) /  canopyi.p_atm
    leaf_gsw1 = leaf_e1 / leaf_d
    leaf_gsc1 = leaf_gsw1 / FTYP(1.6) / ( 1 + canopyi.g_ias_c * leaf_gsw1^(canopyi.g_ias_e) )
    anagrpi   = get_leaf_an_ag_r_pi_from_gsc(
                                             v25 = canopyi.v_max,
                                             j25 = canopyi.j_max,
                                          Γ_star = canopyi.Γ_star,
                                             gsc = leaf_gsc1,
                                             p_a = canopyi.p_a,
                                             tem = canopyi.t_list[indx],
                                             par = canopyi.par_list[indx],
                                           p_atm = canopyi.p_atm,
                                            p_O₂ = canopyi.p_O₂,
                                             r25 = canopyi.r_25)
    leaf_a1   = anagrpi[1]

    # mol CO₂ mol⁻¹ H₂O, Δe = 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return leaf_a1 - canopyi.an_list[indx]
end




"""
    get_marginal_gain(canopyi, indx)

A list of marginal water use efficiency for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
"""
function get_marginal_gain(canopyi::CanopyLayer)
    FTYP      = eltype(canopyi.an_list[1])
    leaf_e1   = canopyi.e_list .+ FTYP(1E-6)
    leaf_gsw1 = leaf_e1 ./ canopyi.d_list
    leaf_gsc1 = leaf_gsw1 ./ FTYP(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.gsw_list .^ (canopyi.g_ias_e) )
    anagrpi_l = get_leaf_an_ag_r_pi_from_gsc_list(
                                                  v25 = canopyi.v_max,
                                                  j25 = canopyi.j_max,
                                               Γ_star = canopyi.Γ_star,
                                             gsc_list = leaf_gsc1,
                                                  p_a = canopyi.p_a,
                                             tem_list = canopyi.t_list,
                                             par_list = canopyi.par_list,
                                                p_atm = canopyi.p_atm,
                                                 p_O₂ = canopyi.p_O₂,
                                                  r25 = canopyi.r_25)
    leaf_a1   = anagrpi_l[1]

    # mol CO₂ mol⁻¹ H₂O, 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return leaf_a1 .- canopyi.an_list
end
