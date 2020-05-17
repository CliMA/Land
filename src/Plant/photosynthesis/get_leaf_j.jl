"""
    get_leaf_j(jmax, par)
This function calculates j from par using a quadratic function
"""
function get_leaf_j(jmax::FT, par::FT) where {FT}
    b = PS_J_QY * par + jmax
    c = PS_J_QY * par * jmax
    j = ( b - sqrt(b^NUMB_2 - NUMB_4*PS_J_CR*c) ) / PS_J_CR * NUMB_0_5
    return j
end
