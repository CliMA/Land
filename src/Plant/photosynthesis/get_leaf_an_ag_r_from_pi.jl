"""
    get_leaf_an_ag_r_from_pi(v25, j25, Γ_star, p_i, tem, par, p_O₂, r25)

Net phtosynthetic rate `an`, gross phtosynthetic rate `ag`, and respiration rate `rday`, given
- `v25` Maixmal carboxylation rate at 298.15 K (25 Celcius)
- `j25` Maximal electron transport rate at 298.15 K (25 Celcius)
- `Γ_star` CO₂ compensation point with the absense of dark respiration
- `p_i` Leaf internal CO₂ partial pressure
- `tem` Leaf temperature
- `par` Photosynthetic active radiation
- `p_O₂` O₂ partial pressure
- `r25` Leaf respiration rate at 298.15 K (25 Celcius). If ==Inf, r25 will be computed from v25.
"""
function get_leaf_an_ag_r_from_pi(;
                                  v25::FT = FT(80.0),
                                  j25::FT = FT(135.0),
                               Γ_star::FT = FT(2.5),
                                  p_i::FT = FT(30.0),
                                  tem::FT = FT(298.15),
                                  par::FT = FT(1000.0),
                                 p_O₂::FT = FT(21278.25),
                                  r25::FT = FT(Inf) ) where {FT}
    # calculate respiration rate
    if r25==Inf
        rday = get_leaf_r_from_v25(v25, tem)
    else
        rday = get_leaf_r_from_r25(r25, tem)
    end

    # calculate ac and aj
    vmax = get_leaf_vcmax(v25, tem)
    jmax = get_leaf_jmax(j25, tem)
    j    = get_leaf_j(jmax,par)
    kc   = PS_KC_25 * PS_KC_Q10^( FT(0.1) * (tem-K_25) )
    ko   = PS_KO_25 * PS_KO_Q10^( FT(0.1) * (tem-K_25) )
    km   = kc * (1 + p_O₂/ko)
    aj   = j * (p_i-Γ_star) / ( 4*(p_i + 2*Γ_star) )
    ac   = vmax * (p_i-Γ_star) / (p_i+km)
    ag   = min(ac,aj)
    an   = ag - rday

    # return a_net
    return an,ag,rday
end
