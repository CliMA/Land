"""
    get_jmax(j25::FT, para_set::AbstractJmax, t_leaf::FT)

Maximal electron transport rate at leaf temperature, given
- `j25` Maximal electron transport rate at 298.15 K
- `para_set` AbstractJmax parameter set
- `t_leaf` Leaf temperature
"""
function get_jmax(j25::FT, para_set::AbstractJmax, t_leaf::FT) where {FT}
    return j25 * arrhenius_correction(para_set, t_leaf)
end
