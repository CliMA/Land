#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function leaf_photosynthesis!
#
#######################################################################################################################################################################################################
"""
This function calculates leaf photosynthetic rates and save the values within the structure. Supported methods are

$(METHODLIST)

"""
function leaf_photosynthesis! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: set a default p_i from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: add examples to docs
#     2022-Jan-14: use colimit function to compute a_gross and a_net
# To do
#     TODO: update leaf T in StomataModels module or higher level
#
#######################################################################################################################################################################################################
"""
    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `p_i` Internal CO₂ partial pressure in `Pa`, default is `leaf.p_CO₂_i`

---
# Examples
```julia
leaf = Leaf{Float64}();
air = AirLayer{FT}();
mode = PCO₂Mode();
leaf_photosynthesis!(leaf, air, mode);
leaf_photosynthesis!(leaf, air, mode, 30.0);
```
"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat} = (
    leaf.p_CO₂_i = p_i;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaf.t != leaf._t
        photosystem_temperature_dependence!(leaf.PSM, air, leaf.t);
    end;
    photosystem_electron_transport!(leaf.PSM, leaf.PRC, leaf.apar);
    rubisco_limited_rate!(leaf.PSM, leaf.p_CO₂_i);
    light_limited_rate!(leaf.PSM, leaf.p_CO₂_i);
    product_limited_rate!(leaf.PSM, leaf.p_CO₂_i);
    colimit_photosynthesis!(leaf.PSM, leaf.COLIMIT);

    # update the fluorescence related parameters
    photosystem_coefficients!(leaf);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: set a default g_lc from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: add examples to docs
#     2022-Jan-14: use colimit function to compute a_gross and a_net
# To do
#     TODO: update leaf T in StomataModels module or higher level
#
#######################################################################################################################################################################################################
"""
    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, g_lc::FT = leaf.g_CO₂) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, default is `leaf.g_CO₂`

---
# Examples
```julia
leaf = Leaf{Float64}();
air = AirLayer{FT}();
mode = GCO₂Mode();
leaf_photosynthesis!(leaf, air, mode);
leaf_photosynthesis!(leaf, air, mode, 0.1);
```
"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, g_lc::FT = leaf.g_CO₂) where {FT<:AbstractFloat} = (
    leaf.g_CO₂ = g_lc;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaf.t != leaf._t
        photosystem_temperature_dependence!(leaf.PSM, air, leaf.t);
    end;
    photosystem_electron_transport!(leaf.PSM, leaf.PRC, leaf.apar);
    rubisco_limited_rate!(leaf.PSM, air, leaf.g_CO₂);
    light_limited_rate!(leaf.PSM, air, leaf.g_CO₂);
    product_limited_rate!(leaf.PSM, air, leaf.g_CO₂);
    colimit_photosynthesis!(leaf.PSM, leaf.COLIMIT);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf.p_CO₂_i = air.p_CO₂ - leaf.PSM.a_net / leaf.g_CO₂   * air.P_AIR * FT(1e-6);
    leaf.p_CO₂_s = air.p_CO₂ - leaf.PSM.a_net / leaf.g_CO₂_b * air.P_AIR * FT(1e-6);

    # update the fluorescence related parameters
    photosystem_coefficients!(leaf);

    return nothing
);
