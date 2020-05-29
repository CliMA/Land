"""
    get_kpep(para_set::AbstractKpep, t_leaf::FT)

Kpep at leaf temperature, given
- `para_set` One `AbstractKpep` type that store temperature correction information
- `t_leaf` Leaf temperature
"""
function get_kpep(paraset::AbstractKpepTD, t_leaf::FT) where {FT}
    return paraset.Kpep * arrhenius_correction(paraset, t_leaf)
end
