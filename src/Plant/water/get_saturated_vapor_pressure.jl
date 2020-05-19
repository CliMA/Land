"""
    get_saturated_vapor_pressure(tem)

Saturated vapor pressure in `[Pa]`, given
- `tem` Water temperature

Saturated vapor pressure is computed using 611.0 * exp(17.502 * temc / (temc + 240.97)).

May need to merge with other CLIMA repository to be consistent.
"""
function get_saturated_vapor_pressure(tem::FT) where {FT}
    temc = tem - FT(K_0)
    return FT(svp_a0) * exp(FT(svp_e0) * temc / (temc + FT(svp_e1)))
end




"""
    get_saturated_vapor_pressure(tem)

List of saturated vapor pressure in `[Pa]`, given
- `tem` A list of temperatures

Saturated vapor pressure is computed using 611.0 * exp(17.502 * temc / (temc + 240.97)).

May need to merge with other CLIMA repository to be consistent.
"""
function get_saturated_vapor_pressure(tem::Array{FT,1}) where {FT}
    temc = tem .- FT(K_0)
    return FT(svp_a0) .* exp.(FT(svp_e0) .* temc ./ (temc .+ FT(svp_e1)))
end
