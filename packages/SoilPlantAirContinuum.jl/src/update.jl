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
#     2022-Oct-19: add method to update or prescribe cab, car, lai, swcs, Vcmax and Jmax TD, t_leaf, vcmax profile
#     2022-Oct-19: air.rh and air.p_H₂O_sat have been removed in an earlier version ClimaCache
#     2022-Oct-19: add method to prescribe t_soil profile
#     2022-Nov-21: fix a bug related to Vcmax profile (no global simulations are impacted)
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
    if !isnothing(t) air.t = t; end;
    if !isnothing(wind) air.wind = wind; end;
    if !isnothing(p_CO₂) air.p_CO₂ = p_CO₂; end;
    if !isnothing(p_H₂O) air.p_H₂O = p_H₂O; end;
    if !isnothing(rh) air.p_H₂O = saturation_vapor_pressure(air.t) * rh; end;
    if !isnothing(vpd) air.p_H₂O = max(0, saturation_vapor_pressure(air.t) - vpd); end;

    return nothing
);


"""

    update!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}};
            cab::Union{Number,Nothing} = nothing,
            car::Union{Number,Nothing} = nothing,
            lai::Union{Number,Nothing} = nothing,
            swcs::Union{Tuple,Nothing} = nothing,
            t_clm::Union{Number,Nothing} = nothing,
            t_leaf::Union{Number,Nothing} = nothing,
            t_soils::Union{Tuple,Nothing} = nothing,
            vcmax::Union{Number,Nothing} = nothing,
            vcmax_expo::Union{Number,Nothing} = nothing
    ) where {FT<:AbstractFloat}

Update the physiological parameters of the SPAC, given
- `spac` Soil plant air continuum
- `cab` Chlorophyll content. Optional, default is nothing
- `car` Carotenoid content. Optional, default is nothing
- `lai` Leaf area index. Optional, default is nothing
- `swcs` Soil water content at different layers. Optional, default is nothing
- `t_clm` Moving average temperature to update Vcmax and Jmax temperature dependencies. Optional, default is nothing
- `t_leaf` Leaf temperature. Optional, default is nothing
- `t_soils` Soil temperature at different layers. Optional, default is nothing
- `vcmax` Vcmax25 at the top of canopy. Optional, default is nothing
-`vcmax_expo` Exponential tuning factor to adjust Vcmax25. Optional, default is nothing

"""
update!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}};
        cab::Union{Number,Nothing} = nothing,
        car::Union{Number,Nothing} = nothing,
        lai::Union{Number,Nothing} = nothing,
        swcs::Union{Tuple,Nothing} = nothing,
        t_clm::Union{Number,Nothing} = nothing,
        t_leaf::Union{Number,Nothing} = nothing,
        t_soils::Union{Tuple,Nothing} = nothing,
        vcmax::Union{Number,Nothing} = nothing,
        vcmax_expo::Union{Number,Nothing} = nothing
) where {FT<:AbstractFloat} = (
    @unpack CANOPY, DIM_LAYER, LEAVES, SOIL = spac;

    # update chlorophyll and carotenoid contents (and )
    if !isnothing(cab)
        for _leaf in LEAVES
            _leaf.BIO.cab = cab;
            _leaf.BIO._v_storage = 0;
        end;
    end;
    if !isnothing(car)
        for _leaf in LEAVES
            _leaf.BIO.car = car;
            _leaf.BIO._v_storage = 0;
        end;
    end;
    if !isnothing(cab) || !isnothing(car)
        leaf_spectra!(spac);
    end;

    # update LAI
    if !isnothing(lai)
        CANOPY.lai = lai;
        for _i in 1:DIM_LAYER
            LEAVES[_i].HS.AREA = SOIL.AREA * CANOPY.lai / DIM_LAYER;
        end;
    end;

    # prescribe soil water content
    if !isnothing(swcs)
        for _i in eachindex(swcs)
            SOIL.LAYERS[_i].θ = swcs[_i];
        end;
    end;

    # prescribe soil temperature
    if !isnothing(swcs)
        for _i in eachindex(swcs)
            SOIL.LAYERS[_i].t = t_soils[_i];
        end;
    end;

    # update Vcmax and Jmax TD
    if !isnothing(t_clm)
        for _leaf in LEAVES
            _leaf.PSM.TD_VCMAX.ΔSV = 668.39 - 1.07 * (t_clm - T₀());
            _leaf.PSM.TD_JMAX.ΔSV = 659.70 - 0.75 * (t_clm - T₀());
        end;
    end;

    # prescribe leaf temperature
    if !isnothing(t_leaf)
        for _leaf in LEAVES
            _leaf.t = t_leaf;
        end;
    end;

    # update Vcmax at the top layer
    if !isnothing(vcmax)
        LEAVES[1].PSM.v_cmax25 = vcmax;
    end;

    # update Vcmax profile if lai or vcmax is given
    if !isnothing(vcmax) || !isnothing(lai)
        for _i in 2:DIM_LAYER
            _scaling = isnothing(vcmax_expo) ? 1 : exp(-vcmax_expo * CANOPY.lai * ((_i - 1) / DIM_LAYER));
            LEAVES[_i].PSM.v_cmax25 = LEAVES[1].PSM.v_cmax25 * _scaling;
            LEAVES[_i].PSM.j_max25 = LEAVES[1].PSM.v_cmax25 * 1.67 * _scaling;
            LEAVES[_i].PSM.r_d25 = LEAVES[1].PSM.v_cmax25 * 0.015 * _scaling;
            LEAVES[_i].PSM._t = 0;
        end;
    end;

    return nothing
);
