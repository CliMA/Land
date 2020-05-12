# calculate j from par
function get_leaf_j(jmax::Number, par::Number)
    a =  0.9
    b = -0.3 * par - jmax
    c =  0.3 * par * jmax
    j = ( -b - sqrt(b*b-4*a*c) ) / a * 0.5
    return j
end
