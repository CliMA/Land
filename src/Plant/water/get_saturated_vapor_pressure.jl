"""
    get_saturated_vapor_pressure(tem)
This function returns the saturated vapor pressure in Pa.
p_sat is computed using 611.0 * exp(17.502 * temc / (temc + 240.97)).
The a0, e0, and e1 are moved to constants
"""
function get_saturated_vapor_pressure(tem::FT) where {FT}
    temc = tem - K_0
    return svp_a0 * exp(svp_e0 * temc / (temc + svp_e1))
end




"""
    get_saturated_vapor_pressure(tem)
This function returns the saturated vapor pressure in Pa for a list of temperatures.
p_sat is computed using 611.0 * exp(17.502 * temc / (temc + 240.97)).
The a0, e0, and e1 are moved to constants
"""
function get_saturated_vapor_pressure(tem::Array{FT,1}) where {FT}
    temc = tem .- K_0
    return svp_a0 .* exp.(svp_e0 .* temc ./ (temc .+ svp_e1))
end
