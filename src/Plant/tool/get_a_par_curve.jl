"""
    get_a_par_curve(v25, j25, Γ_star, gsc, p_a, tem, p_O₂, r25)

Lists of phtosynthetic active radiation `list_par`, net phtosynthetic rate `list_an`, gross phtosynthetic rate `list_ag`, and leaf internal CO₂ partial pressure `list_pi`, given constant
- `v25` Maximal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K
- `Γ_star` CO₂ compensation point with the absence of dark respiration
- `gsc` Leaf diffusive conductance to CO₂
- `p_a` Atmospheric CO₂ partial pressure
- `tem` Leaf temperature
- `p_atm` Atmospheric pressure
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K. If ==Inf, use v25 to compute r25.

# Description
This function is meant to obtain the Anet-PAR curve.
This function may also be used to test the phtosynthesis module.
"""
function get_a_par_curve(model::AbstractPhotoModelParaSet,
                           v25::FT,
                           j25::FT,
                           p25::FT,
                           gsc::FT,
                           p_a::FT,
                           tem::FT,
                         p_atm::FT,
                          p_O₂::FT,
                           r25::FT,
                     curvature::FT,
                            qy::FT) where {FT}
    # generate a list of p_i from gamma to 200 Pa
    list_par = [FT(val) for val in 5.0:5.0:2000.0]
    list_ag  = list_par .* FT(-Inf)
    list_an  = list_par .* FT(-Inf)
    list_pi  = list_par .* FT(-Inf)

    # iterate through the list
    for i in 1:length(list_par)
        par         = list_par[i]
        an,ag,r,p_i = get_an_ag_r_pi_from_gsc(
                                              model,
                                              gsc,
                                              v25,
                                              j25,
                                              p25,
                                              p_a,
                                              tem,
                                              par,
                                              p_atm,
                                              p_O₂,
                                              r25,
                                              curvature,
                                              qy)
        list_ag[i]  = ag
        list_an[i]  = an
        list_pi[i]  = p_i
    end

    # return the arrays
    return list_par,list_an,list_ag,list_pi
end
