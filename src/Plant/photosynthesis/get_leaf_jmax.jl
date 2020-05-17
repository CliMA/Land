"""
    get_leaf_jmax(v25, tem)
This function computes the Jmax from Jmax@25 Celcius
"""
function get_leaf_jmax(jmax25::FT, tem::FT) where {FT}
    factor = PS_J_C * exp( PS_J_HA / GAS_R / K_25 * (NUMB_1 - K_25/tem) ) / ( NUMB_1 + exp( (PS_J_SV*tem-PS_J_HD) / (GAS_R*tem) ) )
    return jmax25 * factor
end
