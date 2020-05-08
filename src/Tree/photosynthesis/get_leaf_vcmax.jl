# function to get vcmax from leaf tempearture
function get_leaf_vcmax(vcmax25, tem, unit="K")
    # convert to degree C
    if unit=="K" || unit=="k"
        temk = tem
    else
        temk = tem + 273.15
    end

    # compute the multiplier
    ha     = 73637.0
    hd     = 149252.0
    sv     = 486.0
    t0     = 298.15
    r      = 8.315
    c      = 1.0 + exp((sv*t0 -hd)/(r*t0))
    factor = c * exp(ha/r/t0*(1.0-t0/temk)) / (1.0 + exp((sv*temk-hd)/(r*temk)))

    # return the vcmax at new temperature
    return vcmax25 * factor
end