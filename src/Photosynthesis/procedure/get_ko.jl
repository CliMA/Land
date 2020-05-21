"""
    get_ko(para_set::AbstractKo, t_leaf::FT)

Ko at leaf temperature, given
- `para_set` One `AbstractKo` type that store temperature correction information
- `t_leaf` Leaf temperature
"""
function get_ko(paraset::AbstractKoTD, t_leaf::FT) where {FT}
    return paraset.Ko * arrhenius_correction(paraset, t_leaf)
end
