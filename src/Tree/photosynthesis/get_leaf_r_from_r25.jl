# get respiration rate from r25 and leaf temperature
function get_leaf_r_from_r25(r25, tem, unit="K")
    # convert to degree C
    if unit=="K" || unit=="k"
        temc = tem - 273.15
    else
        temc = tem
    end

    # compute r from temperature
    r_leaf  = r25 * 2.0^((temc-25.0)*0.1)
    r_leaf /= 1.0 + exp(1.3*(temc-55.0))

    # return respiration rate
    return r_leaf
end