"""
    update_empirical_gsw_empiricla!(tree::Tree, model::AbstractEmpiricalStomatalModel)

Update the gsw target from empirical approach, given
- `tree` A [`Tree`](@ref) type struct
- `model` An empirical stomatal model option.

Current function works for Ball-Berry, Leuning, and Medlyn type models. 


"""
function update_empirical_gsw_empirical!(tree::Tree, model::AbstractEmpiricalStomatalModel)
    # 1. iterate through the canopy layers
    for canopyi in tree.canopy_list
        # 2. iterate through the leaves
        for indx in 1:length(canopyi.leaf_list)
            gsw = get_empirical_gsw(
                                    canopyi.curvature,
                                    canopyi.g_ias_c,
                                    canopyi.g_ias_e,
                                    canopyi.g_max,
                                    canopyi.j_max,
                                    canopyi.p_max,
                                    canopyi.p_a,
                                    canopyi.p_atm,
                                    canopyi.p_H₂O,
                                    canopyi.p_O₂,
                                    canopyi.par_list[indx],
                                    canopyi.qy,
                                    canopyi.r_25,
                                    canopyi.t_list[indx],
                                    canopyi.v_max,
                                    tree.photo_para_set,
                                    model)
            canopyi.gsw_empi[indx] = gsw
        end
    end
end