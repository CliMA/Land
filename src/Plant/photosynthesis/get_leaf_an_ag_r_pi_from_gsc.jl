# get lists of an,ag,r,pi from gsc
function get_leaf_an_ag_r_pi_from_gsc(;
                                       v25::Number = 80.0,
                                       j25::Number = 135.0,
                                    Γ_star::Number = 2.5,
                                       gsc::Number = 0.1,
                                       p_a::Number = 40.0,
                                       tem::Number = 298.15,
                                       par::Number = 1000.0,
                                     p_atm::Number = 101325.0,
                                      p_O₂::Number = 21278.25,
                                       r25::Number = Inf,
                                      unit::String = "K")
    # compute a_net using bi-section method
    tar_p  = 0.0
    tar_ag = 0.0
    tar_an = 0.0
    tar_r  = 0.0
    max_p  = p_a
    min_p  = Γ_star
    while true
        tar_p = 0.5 * (max_p+min_p)
        an,ag,r = get_leaf_an_ag_r_from_pi( v25 = v25,
                                            j25 = j25,
                                         Γ_star = Γ_star,
                                            p_i = tar_p,
                                            tem = tem,
                                            par = par,
                                           p_O₂ = p_O₂,
                                            r25 = r25,
                                           unit = unit)
        tmp_g = an * 1e-6 / (p_a-tar_p) * p_atm

        # increase min_p when g is smaller than target
        if abs(tmp_g-gsc) < 1E-10
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        elseif tmp_g<gsc
            min_p = tar_p
        else
            max_p = tar_p
        end

        # if the max_p and min_p equals, break; used for the case when g=0
        if abs(max_p-min_p) < 1E-6
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        end
    end

    # return the an, ag, r, and p_i
    return tar_an, tar_ag, tar_r, tar_p
end




# get lists of an,ag,r,pi from gsc
function get_leaf_an_ag_r_pi_from_gsc_list(;
                                           v25::Number = 80.0,
                                           j25::Number = 135.0,
                                        Γ_star::Number = 2.5,
                                      gsc_list::Array  = 0.1 .* ones(10),
                                           p_a::Number = 40.0,
                                      tem_list::Array  = 298.15 .* ones(10),
                                      par_list::Array  = 1000.0 .* ones(10),
                                         p_atm::Number = 101325.0,
                                          p_O₂::Number = 21278.25,
                                           r25::Number = Inf,
                                          unit::String = "K")
    # define lists of results
    len     = length(gsc_list)
    list_an = zeros(len)
    list_ag = zeros(len)
    list_re = zeros(len)
    list_pi = zeros(len)

    # iterate the list
    for indx in 1:length( len )
        gsc         = gsc_list[indx]
        tem         = tem_list[indx]
        par         = par_list[indx]
        an,ag,r,p_i = get_leaf_an_ag_r_pi_from_gsc(v25 = v25,
                                                   j25 = j25,
                                                Γ_star = Γ_star,
                                                   gsc = gsc,
                                                   p_a = p_a,
                                                   tem = tem,
                                                   par = par,
                                                 p_atm = p_atm,
                                                  p_O₂ = p_O₂,
                                                   r25 = r25,
                                                  unit = unit)
        list_an[indx] = an
        list_ag[indx] = ag
        list_re[indx]  = r
        list_pi[indx] = p_i
    end

    # return the lists
    return list_an, list_ag, list_re, list_pi
end
