# this function returns the saturated vapor pressure in Pa
function get_saturated_vapor_pressure(tem, unit="K")
    if unit=="c" || unit=="C"
        temc = tem
    else
        temc = tem - 273.15
    end
    
    return 611.0 .* exp.(17.502 .* temc ./ (temc .+ 240.97))
end
