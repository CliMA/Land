#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Apr-19: separate this function as an individual step of the SPAC module (1st step)
#     2022-Jul-12: rename function to update!
#
#######################################################################################################################################################################################################
"""
This function updates the environmental conditions and the soil-plant-air-continuum. Supported functionalities are for
- AirLayer
- SoilPlantAirContinuum

"""
function update! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Apr-19: add the method to update the dignostic variables from air temperature
#     2022-Apr-19: add options p_CO₂, p_H₂O, rh, t, vpd, and wind (defaults are nothing)
#     2022-Apr-19: update docs and history log
#     2022-Jul-12: rename function to update!
#     2022-Jul-12: remove FT control to options
#
#######################################################################################################################################################################################################
"""

    update!(air::AirLayer{FT};
            p_CO₂::Union{Number,Nothing} = nothing,
            p_H₂O::Union{Number,Nothing} = nothing,
            rh::Union{Number,Nothing} = nothing,
            t::Union{Number,Nothing} = nothing,
            vpd::Union{Number,Nothing} = nothing,
            wind::Union{Number,Nothing} = nothing
    ) where {FT<:AbstractFloat}

Update the environmental conditions (such as saturated vapor pressure and relative humidity) of the air surrounding the leaf, given
- `air` `AirLayer` type structure
- `p_CO₂` CO₂ partial pressure in `Pa`. Optional, default is nothing
- `p_H₂O` Vapor pressure in `Pa`. Optional, default is nothing
- `rh` Relatibe humidity (fraction). Optional, default is nothing
- `t` Air temperature in `K`. Optional, default is nothing
- `vpd` Vapor pressure deficit `Pa`. Optional, default is nothing
- `wind` Wind speed in `m s⁻¹`. Optional, default is nothing

"""
update!(air::AirLayer{FT};
        p_CO₂::Union{Number,Nothing} = nothing,
        p_H₂O::Union{Number,Nothing} = nothing,
        rh::Union{Number,Nothing} = nothing,
        t::Union{Number,Nothing} = nothing,
        vpd::Union{Number,Nothing} = nothing,
        wind::Union{Number,Nothing} = nothing
) where {FT<:AbstractFloat} = (
    if !isnothing(t)     air.t = t;         end;
    if !isnothing(p_CO₂) air.p_CO₂ = p_CO₂; end;
    if !isnothing(wind)  air.wind = wind;   end;
    if !isnothing(p_H₂O) air.p_H₂O = p_H₂O; end;

    air.p_H₂O_sat = saturation_vapor_pressure(air.t);

    if !isnothing(rh)  air.p_H₂O = air.p_H₂O_sat * rh;  end;
    if !isnothing(vpd) air.p_H₂O = air.p_H₂O_sat - vpd; end;

    air.rh = air.p_H₂O / air.p_H₂O_sat;

    return nothing
);
