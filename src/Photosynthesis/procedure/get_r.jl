"""
    get_r(para_set::AbstractRespirationTD, r25::FT, t_leaf::FT)

Leaf respiration rate `r_leaf` at leaf temperature, given
- `para_set` One `AbstractRespirationTD` type that store temperature correction information
- `r25` Leaf respiration rate at 298.15 K (25 Celcius)
- `t_leaf` Leaf temperature
"""
function get_r(paraset::AbstractRespirationTD, r25::FT, t_leaf::FT) where {FT}
    return r25 * arrhenius_correction(paraset, t_leaf)
end
