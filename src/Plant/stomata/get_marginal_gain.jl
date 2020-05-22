"""
    get_marginal_gain(canopyi, indx, photo_para_set)

Marginal water use efficiency `∂A/∂E` for a leaf in a canopy layer, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
- `indx` indx'th [`Leaf`](@ref) in the `canopyi.leaf_list`
"""
function get_marginal_gain(canopyi::CanopyLayer, indx::Number; photo_para_set::AbstractPhotoModelParaSet=C3VcVpJBernacchi{FT}()) where {FT}
    leaf_e1   = canopyi.e_list[indx] + FT(1E-6)
    leaf_d    = (get_saturated_vapor_pressure(canopyi.t_list[indx]) - canopyi.p_H₂O) /  canopyi.p_atm
    leaf_gsw1 = leaf_e1 / leaf_d
    leaf_gsc1 = leaf_gsw1 / FT(1.6) / ( 1 + canopyi.g_ias_c * leaf_gsw1^(canopyi.g_ias_e) )
    anagrpi   = get_an_ag_r_pi_from_gsc(
                                        photo_para_set,
                                        leaf_gsc1,
                                        canopyi.v_max,
                                        canopyi.j_max,
                                        canopyi.p_max,
                                        canopyi.p_a,
                                        canopyi.t_list[indx],
                                        canopyi.par_list[indx],
                                        canopyi.p_atm,
                                        canopyi.p_O₂,
                                        canopyi.r_25,
                                        canopyi.curvature,
                                        canopyi.qy)
    leaf_a1   = anagrpi[1]

    # mol CO₂ mol⁻¹ H₂O, Δe = 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return leaf_a1 - canopyi.an_list[indx]
end




"""
    get_marginal_gain(canopyi, indx)

A list of marginal water use efficiency for all the leaves, given
- `canopyi` A [`CanopyLayer`](@ref) in [`Canopy`](@ref) in a [`Tree`](@ref)
"""
function get_marginal_gain(canopyi::CanopyLayer; photo_para_set::AbstractPhotoModelParaSet=C3VcVpJBernacchi{FT}()) where {FT}
    leaf_e1   = canopyi.e_list .+ FT(1E-6)
    leaf_gsw1 = leaf_e1 ./ canopyi.d_list
    leaf_gsc1 = leaf_gsw1 ./ FT(1.6) ./ ( 1 .+ canopyi.g_ias_c .* canopyi.gsw_list .^ (canopyi.g_ias_e) )
    anagrpi_l = get_an_ag_r_pi_from_gsc_list(
                                             photo_para_set,
                                             leaf_gsc1,
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
    leaf_a1   = anagrpi_l[1]

    # mol CO₂ mol⁻¹ H₂O, 1e-6 mol (H₂O) and μmol = 1e-6 mol (CO₂) cancel out
    return leaf_a1 .- canopyi.an_list
end
