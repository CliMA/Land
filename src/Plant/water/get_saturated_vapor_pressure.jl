# this function returns the saturated vapor pressure in Pa
function get_saturated_vapor_pressure(tem::Number; unit::String="K")
    if unit=="c" || unit=="C"
        temc = tem
    else
        temc = tem - 273.15
    end
    
    return 611.0 * exp.(17.502 * temc / (temc + 240.97))
end




# this function returns the saturated vapor pressure in Pa for a list of temperatures
function get_saturated_vapor_pressure(tem::Array; unit::String="K")
    if unit=="c" || unit=="C"
        temc = tem
    else
        temc = tem .- 273.15
    end
    
    return 611.0 .* exp.(17.502 .* temc ./ (temc .+ 240.97))
end
