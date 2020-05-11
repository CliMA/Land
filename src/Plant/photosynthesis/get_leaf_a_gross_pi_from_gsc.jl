# get ci and Anet from gc and ca
function get_leaf_a_gross_pi_from_gsc(v25, j25, gamma, gsc, p_a, tem, par, p_atm=101325.0, p_o2=21278.25, r25=false, unit="K")
    # calculate the respiration rate
    if r25==false
        rday = get_leaf_r_from_v25(v25, tem, unit)
    else
        rday = get_leaf_r_from_r25(r25, tem, unit)
    end

    # compute a_net using bi-section method
    tar_a,tar_p = get_leaf_a_net_pi_from_gsc(v25, j25, gamma, gsc, p_a, tem, par, p_atm, p_o2, r25, unit)

    # return the a and p
    return tar_a+rday, tar_p
end
