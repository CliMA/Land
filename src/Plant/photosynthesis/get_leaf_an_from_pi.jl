# get a from ci
function get_leaf_an_from_pi(v25, j25, gamma, p_i, tem, par, p_o2=21278.25, r25=false, unit="K")
    # convert to degree C
    if unit=="K" || unit=="k"
        temk = tem
    else
        temk = tem + 273.15
    end

    # calculate respiration rate
    if r25==false
        rday = get_leaf_r_from_v25(v25, temk)
    else
        rday = get_leaf_r_from_r25(r25, temk)
    end

    # calculate ac and aj
    vmax = get_leaf_vcmax(v25, temk)
    jmax = get_leaf_jmax(j25, temk)
    j    = get_leaf_j(jmax,par)
    kc   = 41.01637 * 2.1^(0.1*(temk-298.15))
    ko   = 28201.92 * 1.2^(0.1*(temk-298.15))
    km   = kc * (1.0+p_o2/ko)
    aj   = j * (p_i-gamma) / (4.0*(p_i+2*gamma))
    ac   = vmax * (p_i-gamma) / (p_i+km)
    af   = min(ac,aj) - rday

    # return a_net
    return af
end
