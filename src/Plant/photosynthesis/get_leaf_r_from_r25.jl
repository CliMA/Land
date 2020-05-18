"""
    get_leaf_r_from_v25(v25, tem)

# Arguments
- `r25::FT`    Leaf respiration rate at 298.15 K (25 Celcius)
- `tem::FT`    Leaf temperature

# Description
This function returns respiration rate from r25 and leaf temperature
"""
function get_leaf_r_from_r25(r25::FT, tem::FT) where {FT}
    # compute r from temperature in K
    r_leaf  = r25 * PS_R_Q10_E^( (tem-K_25)* FT(0.1) )
    r_leaf /= 1 + exp( PS_R_EXP_K * (tem - PS_R_EXP_T) )

    # return respiration rate
    return r_leaf
end
