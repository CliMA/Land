"""
    get_relative_viscosity(tem)

Viscosity relative to 25 degree C (298.15 K), given
- `tem` Water temperature

Equations used are
υ     = A * exp( B/tem + C*tem + D*tem^2 )
υ/υ25 = exp( B/tem + C*tem + D*tem^2 - B/tem25 - C*tem25 - D*tem25^2 )
fitting parameters are from Reid, Prausnitz, & Poling (1987), valid through 273-643 K
A = 1.856E-11 mPa s
B = 4209      K
C = 0.04527   K⁻¹
D = -3.376E-5 K⁻²

May need to merge with other CLIMA repository to be consistent.
"""
function get_relative_viscosity(tem::FT) where {FT}
    return exp( vis_b * ( 1/tem - 1/K_25) + vis_c * (tem - K_25) + vis_d * (tem^2 - K_25^2) )
end
