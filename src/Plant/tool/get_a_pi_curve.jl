"""
    get_a_pi_curve(v25, j25, Γ_star, tem, par, p_O₂, r25)

Lists of leaf internal CO₂ partial pressure `list_pi`, net photosynthetic rate `list_an`, and gross photosynthetic rate `list_ag`, given constant
- `v25` Maximal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K
- `Γ_star` CO₂ compensation point with the absence of dark respiration
- `tem` Leaf temperature
- `par` Photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K. If ==Inf, use v25 to compute r25.
"""
function get_a_pi_curve(model::C3ParaSet,
                          v25::FT,
                          j25::FT,
                          p25::FT,
                          tem::FT,
                          par::FT,
                         p_O₂::FT,
                          r25::FT,
                    curvature::FT,
                           qy::FT) where {FT}
    # generate a list of p_i from gamma to 200 Pa
    list_pi = [FT(val) for val in 1.0:1.0:200.0]
    list_ag = list_pi .* FT(-Inf)
    list_an = list_pi .* FT(-Inf)

    # iterate through the list
    for i in 1:length(list_pi)
        p_i        = list_pi[i]
        a_n,a_g,r  = get_an_ag_r_from_pi(
                                         model,
                                         p_i,
                                         v25,
                                         j25,
                                         par,
                                         p_O₂,
                                         tem,
                                         r25,
                                         curvature,
                                         qy)
        list_ag[i] = a_g
        list_an[i] = a_n
    end

    # return the arrays
    return list_pi,list_an,list_ag
end
