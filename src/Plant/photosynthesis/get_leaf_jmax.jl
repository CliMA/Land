# function yo get jmax from leaf temperature
function get_leaf_jmax(jmax25::Number, tem::Number; unit::String="K")
    # convert to degree C
    if unit=="K" || unit=="k"
        temk = tem
    else
        temk = tem + 273.15
    end

    # compute the multiplier
    ha     = 50300.0
    hd     = 152044.0
    sv     = 495.0
    t0     = 298.15
    r      = 8.315
    c      = 1.0 + exp((sv*t0 -hd)/(r*t0))
    factor = c * exp(ha/r/t0*(1.0-t0/temk)) / (1.0 + exp((sv*temk-hd)/(r*temk)))

    # return the jmax at new temperature
    return jmax25 * factor
end
