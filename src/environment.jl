#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Apr-19: separate this function as an individual step of the SPAC module (1st step)
#
#######################################################################################################################################################################################################
"""
This function updates the environmental conditions of the soil-plant-air-continuum. Supported methods are

$(METHODLIST)

"""
function update_environment! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Apr-19: add the method to update the dignostic variables from air temperature
#     2022-Apr-19: add options p_H₂O, rh, and t (defaults are nothing)
#
#######################################################################################################################################################################################################
"""

    update_environment!(
                envir::AirLayer{FT};
                p_CO₂::Union{FT,Nothing} = nothing,
                p_H₂O::Union{FT,Nothing} = nothing,
                rh::Union{FT,Nothing} = nothing,
                t::Union{FT, Nothing} = nothing,
                vpd::Union{FT, Nothing} = nothing,
                wind::Union{FT,Nothing} = nothing
    ) where {FT<:AbstractFloat}

Update the environmental conditions (such as saturated vapor pressure and relative humidity) of the air surrounding the leaf based on air temperature, given
- `envir` `AirLayer` type structure
- `p_CO₂` CO₂ partial pressure in `Pa`. Optional, default is nothing
- `p_H₂O` Vapor pressure in `Pa`. Optional, default is nothing
- `rh` Relatibe humidity (fraction). Optional, default is nothing
- `t` Air temperature in `K`. Optional, default is nothing
- `vpd` Vapor pressure deficit `Pa`. Optional, default is nothing
- `wind` Wind speed in `m s⁻¹`. Optional, default is nothing
"""
update_environment!(
            envir::AirLayer{FT};
            p_CO₂::Union{FT,Nothing} = nothing,
            p_H₂O::Union{FT,Nothing} = nothing,
            rh::Union{FT,Nothing} = nothing,
            t::Union{FT, Nothing} = nothing,
            vpd::Union{FT, Nothing} = nothing,
            wind::Union{FT,Nothing} = nothing
) where {FT<:AbstractFloat} = (
    if !isnothing(t)     envir.t = t;         end;
    if !isnothing(p_CO₂) envir.p_CO₂ = p_CO₂; end;
    if !isnothing(wind)  envir.wind = wind;   end;
    if !isnothing(p_H₂O) envir.p_H₂O = p_H₂O; end;

    envir.p_H₂O_sat = saturation_vapor_pressure(envir.t);

    if !isnothing(rh)  envir.p_H₂O = envir.p_H₂O_sat * rh;  end;
    if !isnothing(vpd) envir.p_H₂O = envir.p_H₂O_sat - vpd; end;

    envir.rh = envir.p_H₂O / envir.p_H₂O_sat;

    return nothing
);
