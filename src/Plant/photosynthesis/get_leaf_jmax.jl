"""
    get_leaf_jmax(v25, tem)

# Arguments
- `j25::FT`    Maximal electron transport rate at 298.15 K
- `tem::FT`    Leaf temperature

# DEscription
This function computes the Jmax from Jmax @ 298.15 K (25 Celcius)
"""
function get_leaf_jmax(j25::FT, tem::FT) where {FT}
    factor = PS_J_C * exp( PS_J_HA / GAS_R / K_25 * (1 - K_25/tem) ) / ( 1 + exp( (PS_J_SV*tem-PS_J_HD) / (GAS_R*tem) ) )
    return j25 * factor
end
