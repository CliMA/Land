# this function returns the saturated vapor pressure in Pa
function get_saturated_vapor_pressure(tem::FT) where {FT}
    temc = tem - 273.15
    return 611.0 * exp.(17.502 * temc / (temc + 240.97))
end




# this function returns the saturated vapor pressure in Pa for a list of temperatures
function get_saturated_vapor_pressure(tem::Array)
    temc = tem .- 273.15
    return 611.0 .* exp.(17.502 .* temc ./ (temc .+ 240.97))
end
