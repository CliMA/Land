"""
    get_leaf_j(jmax, par)

Electron transport rate `j`, given
- `jmax` Maximal eclectron transport @ leaf temperature (not 298.15 K)
- `par` Photosynthetic active radiation
"""
function get_leaf_j(jmax::FT, par::FT) where {FT}
    a = FT(PS_J_CR)
    b = FT(PS_J_QY) * par + jmax
    c = FT(PS_J_QY) * par * jmax
    j = ( b - sqrt(b^2 - 4*a*c) ) / (2*a)
    return j
end
