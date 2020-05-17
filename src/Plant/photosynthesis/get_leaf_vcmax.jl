"""
    get_leaf_vcmax(v25, tem)
This function computes the Vcmax from Vcmax@25 Celcius
"""
function get_leaf_vcmax(v25::FT, tem::FT) where {FT}
    factor = PS_V_C * exp( PS_V_HA / GAS_R / K_25 * (NUMB_1 - K_25/tem) ) / ( NUMB_1 + exp( (PS_V_SV*tem-PS_V_HD) / (GAS_R*tem) ) )
    return v25 * factor
end
