"""
    get_vmax(v25::FT, para_set::AbstractVmax, t_leaf::FT)

Maximal electron transport rate at leaf temperature, given
- `para_set` AbstractVmax parameter set
- `v25` Maximal carboxylation rate at 298.15 K
- `t_leaf` Leaf temperature
"""
function get_vmax(para_set::AbstractVmaxTD, v25::FT, t_leaf::FT) where {FT}
    return v25 * arrhenius_correction(para_set, t_leaf)
end
