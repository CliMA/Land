# function yo get jmax from leaf temperature
function get_leaf_jmax(jmax25::FT, tem::FT) where {FT}
    # compute the multiplier
    ha     = 50300.0
    hd     = 152044.0
    sv     = 495.0
    t0     = 298.15
    r      = 8.315
    c      = 1.0 + exp((sv*t0 -hd)/(r*t0))
    factor = c * exp(ha/r/t0*(1.0-t0/tem)) / (1.0 + exp((sv*tem-hd)/(r*tem)))

    # return the jmax at new temperature
    return jmax25 * factor
end
