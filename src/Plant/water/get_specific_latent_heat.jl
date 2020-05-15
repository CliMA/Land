# this function returns the specific latent heat of water
function get_specific_latent_heat(tem::FT) where {FT}
    #=
    lambda = (2500.8 - 2.36*tem + 0.0016*tem^2 -0.00006*tem^3) when tem in -25 to 40 degree C
    see Polynomial curve fits to Table 2.1. R. R. Rogers; M. K. Yau (1989). A Short Course in Cloud Physics (3rd ed.). Pergamon Press. p. 16. ISBN 0-7506-3215-1.
    =#
    temc = tem - 273.15
    return 2500.8 - 2.36*temc + 0.0016*temc^2.0 - 0.00006*temc^3.0
end
