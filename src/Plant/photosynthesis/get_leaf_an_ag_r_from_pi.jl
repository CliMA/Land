# get a from ci
function get_leaf_an_ag_r_from_pi(;
                                  v25::Number = 80.0,
                                  j25::Number = 135.0,
                               Γ_star::Number = 2.5,
                                  p_i::Number = 30.0,
                                  tem::Number = 298.15,
                                  par::Number = 1000.0,
                                 p_O₂::Number = 21278.25,
                                  r25::Number = Inf,
                                 unit::String = "K")
    # convert to degree C
    if unit=="K" || unit=="k"
        temk = tem
    else
        temk = tem + 273.15
    end

    # calculate respiration rate
    if r25==Inf
        rday = get_leaf_r_from_v25(v25, tem; unit=unit)
    else
        rday = get_leaf_r_from_r25(r25, tem; unit=unit)
    end

    # calculate ac and aj
    vmax = get_leaf_vcmax(v25, temk)
    jmax = get_leaf_jmax(j25, temk)
    j    = get_leaf_j(jmax,par)
    kc   = 41.01637 * 2.1^(0.1*(temk-298.15))
    ko   = 28201.92 * 1.2^(0.1*(temk-298.15))
    km   = kc * (1.0+p_O₂/ko)
    aj   = j * (p_i-Γ_star) / (4.0*(p_i+2*Γ_star))
    ac   = vmax * (p_i-Γ_star) / (p_i+km)
    ag   = min(ac,aj)
    an   = ag - rday

    # return a_net
    return an,ag,rday
end
