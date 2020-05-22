"""
    get_j(jmax::FT, par::FT; curvature::FT, f_etr::FT)

Electron transport rate `j`, given
- `jmax` Maximal eclectron transport @ leaf temperature (not 298.15 K)
- `par` Absorbed PAR (photosynthetic active radiation)
- `curvature` Curvature parameter (default at 0.9)
- `qy` Quantum yield of electron, `qy = maxPSII * PSII_frac`, maxPSII is maximal PSII yield (default at 4.0/4.9), PSII_frac is Fraction of absorbed light used by PSII ETR (default at 0.5)
"""
function get_j(jmax::FT, par::FT; curvature::FT=FT(0.9), qy::FT=FT(0.4081632653061224)) where {FT}
    # a = curvature
    b = qy * par + jmax
    c = qy * par * jmax
    j = ( b - sqrt(b^2 - 4*curvature*c) ) / (2*curvature)
    return j
end
