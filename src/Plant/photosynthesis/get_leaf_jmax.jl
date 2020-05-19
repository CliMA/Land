"""
    get_leaf_jmax(v25, tem)

Maximal electron transport rate at leaf temperature, given
- `j25` Maximal electron transport rate at 298.15 K
- `tem` Leaf temperature
"""
function get_leaf_jmax(j25::FT, tem::FT) where {FT}
    factor = FT(PS_J_C) * exp( FT(PS_J_HA) / FT(GAS_R * K_25) * (1 - FT(K_25)/tem) ) / ( 1 + exp( (FT(PS_J_SV)*tem-FT(PS_J_HD)) / (FT(GAS_R)*tem) ) )
    return j25 * factor
end
