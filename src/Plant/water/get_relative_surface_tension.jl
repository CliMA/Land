# this function returns the surface tension of water relative to 25 degree C (298.15 K)
function get_relative_surface_tension(tem::Number; unit::String="K")
    #=
    gamma        = V^(-2/3) * k * (tem_c-tem)
    gamma/gamm25 = (tem_c - tem) / (tem_c - 298.15)
    the empirical values are
    k     = 2.1E-7 J K^-1 mol^(-2/3)
    V     = 18.0 ml/mol
    tem_c = 647.0 K
    =#
    temm = 647.0
    if unit=="C" || unit=="c"
        temk = tem + 273.15
    else
        temk = tem
    end

    return (temm - temk) / (647.0 - 298.15)
end
