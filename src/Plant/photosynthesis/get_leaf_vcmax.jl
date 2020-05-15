# function to get vcmax from leaf tempearture
function get_leaf_vcmax(vcmax25::FT, tem::FT) where {FT}
    # compute the multiplier
    ha     = 73637.0
    hd     = 149252.0
    sv     = 486.0
    t0     = 298.15
    r      = 8.315
    c      = 1.0 + exp((sv*t0 -hd)/(r*t0))
    factor = c * exp(ha/r/t0*(1.0-t0/tem)) / (1.0 + exp((sv*tem-hd)/(r*tem)))

    # return the vcmax at new temperature
    return vcmax25 * factor
end
