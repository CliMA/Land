# this function returns the viscosity relative to 25 degree C (298.15 K)
function get_relative_viscosity(tem::FT) where {FT}
    #=
    mu      = A * exp( B/tem + C*tem + D*tem^2 )
    mu/mu25 = exp( B/tem + C*tem + D*tem^2 - B/tem25 - C*tem25 - D*tem25^2 )
    fitting parameters are from Reid, Prausnitz, & Poling (1987), valid through 273-643 K
    A = 1.856E-11 mPa s
    B = 4209      K
    C = 0.04527   K^-1
    D = -3.376E-5 K^-2
    =#
    B =  4209.0
    C =  0.04527
    D = -3.376E-5

    return exp( B * (1.0/tem - 1.0/298.15) + C * (tem - 298.15) + D * (tem^2.0 - 298.15^2.0) )
end
