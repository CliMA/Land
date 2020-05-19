"""
    get_leaf_r_from_v25(v25, tem)

Leaf respiration rate `r_leaf` at leaf temperature, given
- `r25` Leaf respiration rate at 298.15 K (25 Celcius)
- `tem` Leaf temperature
"""
function get_leaf_r_from_r25(r25::FT, tem::FT) where {FT}
    # compute r from temperature in K
    r_leaf  = r25 * FT(PS_R_Q10_E)^( (tem-FT(K_25))* FT(0.1) )
    r_leaf /= 1 + exp( FT(PS_R_EXP_K) * (tem - FT(PS_R_EXP_T)) )

    # return respiration rate
    return r_leaf
end
