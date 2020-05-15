# get respiration rate from r25 and leaf temperature
function get_leaf_r_from_r25(r25::Number, tem::Number)
    # compute r from temperature in K
    r_leaf  = r25 * 2.0^((tem-298.15)*0.1)
    r_leaf /= 1.0 + exp(1.3*(tem-328.15))

    # return respiration rate
    return r_leaf
end
