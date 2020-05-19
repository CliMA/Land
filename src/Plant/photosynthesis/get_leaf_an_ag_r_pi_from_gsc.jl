"""
    get_leaf_an_ag_r_pi_from_gsc(v25, j25, Γ_star, gsc, p_a, tem, par, p_atm, p_O₂, r25)

Net photosynthetic rate `tar_an`, gross photosynthetic rate `tar_ag`, respiration rate `tar_r`, and leaf internal CO₂ partial pressure `tar_p`, given
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `Γ_star` CO₂ compensation point with the absense of dark respiration
- `gsc` Leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `tem` Leaf temperature
- `par` Photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25.
"""
function get_leaf_an_ag_r_pi_from_gsc(;
                                       v25::FT = FT(80.0),
                                       j25::FT = FT(135.0),
                                    Γ_star::FT = FT(2.5),
                                       gsc::FT = FT(0.1),
                                       p_a::FT = FT(40.0),
                                       tem::FT = FT(298.15),
                                       par::FT = FT(1000.0),
                                     p_atm::FT = FT(101325.0),
                                      p_O₂::FT = FT(21278.25),
                                       r25::FT = FT(Inf) ) where {FT}
    # compute a_net using bi-section method
    tar_p  = FT(0.0)
    tar_ag = FT(0.0)
    tar_an = FT(0.0)
    tar_r  = FT(0.0)
    max_p  = p_a
    min_p  = Γ_star
    while true
        tar_p = FT(0.5) * (max_p+min_p)
        an,ag,r = get_leaf_an_ag_r_from_pi(
                                           v25 = v25,
                                           j25 = j25,
                                        Γ_star = Γ_star,
                                           p_i = tar_p,
                                           tem = tem,
                                           par = par,
                                          p_O₂ = p_O₂,
                                           r25 = r25)
        tmp_g = an * FT(1e-6) / (p_a-tar_p) * p_atm

        # increase min_p when g is smaller than target
        if abs(tmp_g-gsc) < 1e-5
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
        if abs(max_p-min_p) < 1e-5
            tar_ag = ag
            tar_an = an
            tar_r  = r
            break
        end
    end

    # return the an, ag, r, and p_i
    return tar_an, tar_ag, tar_r, tar_p
end




"""
    get_leaf_an_ag_r_pi_from_gsc(v25, j25, Γ_star, gsc_list, p_a, tem_list, par_list, p_atm, p_O₂, r25)

Lists of net photosynthetic rate `list_an`, gross photosynthetic rate `list_ag`, respiration rate `list_re`, and leaf internal CO₂ partial pressure `list_pi`, given
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `Γ_star` CO₂ compensation point with the absense of dark respiration
- `gsc_list` A list of leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `tem_list` A list of leaf temperature
- `par_list` A list of photosynthetic active radiation
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25.
"""
function get_leaf_an_ag_r_pi_from_gsc_list(;
                                           v25::FT          = FT(80.0),
                                           j25::FT          = FT(135.0),
                                        Γ_star::FT          = FT(2.5),
                                      gsc_list::Array{FT,1} = FT(0.1) .* ones(FT,10),
                                           p_a::FT          = FT(40.0),
                                      tem_list::Array{FT,1} = FT(298.15) .* ones(FT,10),
                                      par_list::Array{FT,1} = FT(1000.0) .* ones(FT,10),
                                         p_atm::FT          = FT(101325.0),
                                          p_O₂::FT          = FT(21278.25),
                                           r25::FT          = FT(Inf) ) where {FT}
    # define lists of results
    len     = length(gsc_list)
    list_an = zeros(FT,len)
    list_ag = zeros(FT,len)
    list_re = zeros(FT,len)
    list_pi = zeros(FT,len)

    # iterate the list
    for indx in 1:length( len )
        gsc         = gsc_list[indx]
        tem         = tem_list[indx]
        par         = par_list[indx]
        an,ag,r,p_i = get_leaf_an_ag_r_pi_from_gsc(
                                                   v25 = v25,
                                                   j25 = j25,
                                                Γ_star = Γ_star,
                                                   gsc = gsc,
                                                   p_a = p_a,
                                                   tem = tem,
                                                   par = par,
                                                 p_atm = p_atm,
                                                  p_O₂ = p_O₂,
                                                   r25 = r25)
        list_an[indx] = an
        list_ag[indx] = ag
        list_re[indx] = r
        list_pi[indx] = p_i
    end

    # return the lists
    return list_an, list_ag, list_re, list_pi
end
