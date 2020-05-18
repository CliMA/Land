"""
    get_leaf_r_from_v25(v25, tem)

Leaf respiration rate at leaf temperature if r25 is unknown, given
- `v25` Maximal carboxylation rate at 298.15 K (25 Celcius)
- `tem` Leaf temperature
"""
function get_leaf_r_from_v25(v25::FT, tem::FT) where {FT}
    # compute r25 from v25
    r25 = v25 * PS_RV_RATIO

    # return the r from r25 and temperature
    return get_leaf_r_from_r25(r25, tem)
end
