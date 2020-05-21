"""
    get_kc(para_set::AbstractKc, t_leaf::FT)

Kc at leaf temperature, given
- `para_set` One `AbstractKc` type that store temperature correction information
- `t_leaf` Leaf temperature
"""
function get_kc(paraset::AbstractKcTD, t_leaf::FT) where {FT}
    return paraset.Kc * arrhenius_correction(paraset, t_leaf)
end
