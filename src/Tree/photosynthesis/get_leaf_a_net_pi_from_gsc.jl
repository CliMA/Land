# get ci and Anet from gc and ca
function get_leaf_a_net_ci_from_gsc(v25, j25, gamma, gsc, p_a, tem, par, p_atm=101325.0, p_o2=21000.0, r25=false, unit="K")
    # compute a_net using bi-section method
    tar_p = 0.0
    tar_a = 0.0
    max_p = p_a
    min_p = gamma
    while true
        tar_p = 0.5 * (max_p+min_p)
        af    = get_leaf_a_net_from_pi(v25, j25, gamma, tar_p, tem, par, p_o2, r25, unit)
        tmp_g = af * 1e-6 / (p_a-tar_p) * p_atm

        # increase min_p when g is smaller than target
        if abs(tmp_g-gsc) < 1E-10
            tar_a = af
            break
        elseif tmp_g<gsc
            min_p = tar_p
        else
            max_p = tar_p
        end

        # if the max_p and min_p equals, break; used for the case when g=0
        if abs(max_p-min_p) < 1E-6
            tar_a = af
            break
        end
    end

    # return the a and p
    return tar_a, tar_p
end