# this function returns the saturated vapor pressure in kPa
function get_saturated_vapor_pressure(tem, unit="K")
    if unit=="c" || unit=="C"
        temc = tem
    else
        temc = tem - 273.15
    end
    return 0.611 .* exp.(17.502 .* temc ./ (temc .+ 240.97))
end
