"""
    update_leaf_a_max!(tree::Tree)

Update leaf maximal `A_net` allowed at the given environmental conditions, given
-`tree` A [`Tree`](@ref) type struct

This function is meant to work with Sperry et al. (2017) model, which requires a_max in the ∂Θ∂E function.

"""
function update_leaf_a_max!(tree::Tree)
    @unpack canopy_list,photo_para_set = tree
    FTYP = eltype(tree.ba)

    # 1. iterate through canopy layers
    for canopyi in canopy_list
        # 2. iterate through the leaves in the canopy layers
        for i in 1:length(canopyi.leaf_list)
            leaf = canopyi[i]
            
            # compute a_max
            e_c  = canopyi.ec_list[i]
            gsw  = min( e_c / (canopyi.d_list[i] / canopyi.p_atm), canopyi.g_max )
            gsc  = gsw / FTYP(1.6) / ( 1 + canopyi.g_ias_c * gsw^(canopyi.g_ias_e) )
            re_c = get_an_ag_r_pi_from_gsc(
                                           photo_para_set,
                                           gsc,
                                           canopyi.v_max,
                                           canopyi.j_max,
                                           canopyi.p_max,
                                           canopyi.p_a,
                                           canopyi.t_list[i],
                                           canopyi.par_list[i],
                                           canopyi.p_atm,
                                           canopyi.p_O₂,
                                           canopyi.r_25,
                                           canopyi.curvature,
                                           canopyi.qy)
            a_m  = re_c[1]

            # compute k_max
            p_m  = get_struct_p_end_from_q(leaf, FTYP(0.0))
            k_m  = exp( -(-p_m/leaf.b)^leaf.c )

            # update a_max and k_max
            canopyi.am_list[i] = a_m
            canopyi.km_list[i] = k_m
        end
    end
end