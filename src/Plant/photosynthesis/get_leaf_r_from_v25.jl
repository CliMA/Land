# get respiration rate from vcmax and leaf temperature
function get_leaf_r_from_v25(v25::Number, tem::Number; unit::String="K")
    # compute r25 from v25
    r25 = v25 * 0.015

    # return the r from r25 and temperature
    return get_leaf_r_from_r25(r25, tem; unit=unit)
end
