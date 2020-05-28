"""
    get_jmax(para_set::AbstractJmaxTD, j25::FT, t_leaf::FT)

Maximal electron transport rate at leaf temperature, given
- `para_set` AbstractJmax parameter set
- `j25` Maximal electron transport rate at 298.15 K
- `t_leaf` Leaf temperature
"""
function get_jmax(para_set::AbstractJmaxTD, j25::FT, t_leaf::FT) where {FT}
    return j25 * arrhenius_correction(para_set, t_leaf)
end
