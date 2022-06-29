#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function leaf_photosynthesis!
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""
Per refactored Photosynthesis module, the only things one need to know is the public function `leaf_photosynthesis!` and some construtors from `ClimaCache`. See the examples in the methods below for
    details about how to use the function. The steps for computing photosynthetic rates are

- Update temperature dependent variables using [`photosystem_temperature_dependence!`](@ref)
- Calculate electron transport rate using [`photosystem_electron_transport!`](@ref)
- Calculate RubisCO limited rate using [`rubisco_limited_rate!`](@ref)
- Calculate light limited rate using [`light_limited_rate!`](@ref)
- Calculate product limited rate using [`product_limited_rate!`](@ref)
- Calculate gross and net rates using [`colimit_photosynthesis!`](@ref)
- Update fluorescence related variables using [`photosystem_coefficients!`](@ref)

Supported methods are

$(METHODLIST)

"""
function leaf_photosynthesis! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: set a default p_i from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: add examples to docs
#     2022-Jan-14: use colimit function to compute a_gross and a_net
#     2022-Jan-18: add p_i to electron transport function input variables
#     2022-Jan-24: fix documentation
#     2022-Feb-07: use new method of photosystem_coefficients!
#     2022-Feb-28: use updated light_limited_rate! function
#     2022-Jun-27: use ClimaCache v0.4, where Leaf.apar is renamed to Leaf.ppar
#     2022-Jun-28: unpack the constant fields
# Bug fixes
#     2022-Jan-24: fix PSM abstraction in colimit_photosynthesis! function
# To do
#     TODO: update leaf T in StomataModels module or higher level
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ partial pressure, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `p_i` Internal CO₂ partial pressure in `Pa`, default is `leaf.p_CO₂_i`

---
# Examples
```julia
leaf = Leaf{Float64}("C3");
air  = AirLayer{Float64}();
mode = PCO₂Mode();
leaf_photosynthesis!(leaf, air, mode);
leaf_photosynthesis!(leaf, air, mode, 30.0);
```
"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat} = (
    @unpack PRC, PSM = leaf;

    leaf.p_CO₂_i = p_i;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaf.t != leaf._t
        photosystem_temperature_dependence!(PSM, air, leaf.t);
    end;
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, p_i);
    rubisco_limited_rate!(PSM, leaf.p_CO₂_i);
    light_limited_rate!(PSM);
    product_limited_rate!(PSM, leaf.p_CO₂_i);
    colimit_photosynthesis!(PSM);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaf.ppar);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-28: add method for Leaves1D
#     2022-Jun-28: fix documentation
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::PCO₂Mode) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ partial pressure for sunlit and shaded leaves, given
- `leaves` `Leaves1D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates

---
# Examples
```julia
leaves = Leaves1D{Float64}("C3");
air    = AirLayer{Float64}();
mode   = PCO₂Mode();
leaf_photosynthesis!(leaves, air, mode);
```
"""
leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::PCO₂Mode) where {FT<:AbstractFloat} = (
    @unpack PRC, PSM = leaves;

    # loop through the ppars
    for _i in eachindex(leaves.ppar)
        # update TD parameters everytime for sunlit and shaded leaves
        photosystem_temperature_dependence!(PSM, air, leaves.t[_i]);

        # calculate the photosynthetic rates
        photosystem_electron_transport!(PSM, PRC, leaves.ppar[_i], leaves.p_CO₂_i[_i]);
        rubisco_limited_rate!(PSM, leaves.p_CO₂_i[_i]);
        light_limited_rate!(PSM);
        product_limited_rate!(PSM, leaves.p_CO₂_i[_i]);
        colimit_photosynthesis!(PSM);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar[_i]);

        # save the rates and to leaves
        leaves.a_net[_i] = PSM.a_net;
        leaves.a_gross[_i] = PSM.a_gross;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-28: add method for Leaves1D
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::PCO₂Mode) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ partial pressure for sunlit and shaded leaves, given
- `leaves` `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates

---
# Examples
```julia
leaves = Leaves2D{Float64}("C3");
air    = AirLayer{Float64}();
mode   = PCO₂Mode();
leaf_photosynthesis!(leaves, air, mode);
```
"""
leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::PCO₂Mode) where {FT<:AbstractFloat} = (
    @unpack PRC, PSM = leaves;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaves.t != leaves._t
        photosystem_temperature_dependence!(PSM, air, leaves.t);
    end;

    # loop through the ppars for sunlit leaves
    for _i in eachindex(leaves.ppar_sunlit)
        # calculate the photosynthetic rates
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves.p_CO₂_i_sunlit[_i]);
        rubisco_limited_rate!(PSM, leaves.p_CO₂_i_sunlit[_i]);
        light_limited_rate!(PSM);
        product_limited_rate!(PSM, leaves.p_CO₂_i_sunlit[_i]);
        colimit_photosynthesis!(PSM);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar_sunlit[_i]);

        # save the rates and to leaves
        leaves.a_net_sunlit[_i] = PSM.a_net;
        leaves.a_gross_sunlit[_i] = PSM.a_gross;
        leaves.ϕ_f_sunlit[_i] = PRC.ϕ_f;
    end;

    # run the model for shaded leaves
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves.p_CO₂_i_shaded);
    rubisco_limited_rate!(PSM, leaves.p_CO₂_i_shaded);
    light_limited_rate!(PSM);
    product_limited_rate!(PSM, leaves.p_CO₂_i_shaded);
    colimit_photosynthesis!(PSM);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaves.ppar_shaded);

    # save the rates and to leaves
    leaves.a_net_shaded = PSM.a_net;
    leaves.a_gross_shaded = PSM.a_gross;
    leaves.ϕ_f_shaded = PRC.ϕ_f;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: set a default g_lc from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: add examples to docs
#     2022-Jan-14: use colimit function to compute a_gross and a_net
#     2022-Jan-24: fix documentation
#     2022-Feb-07: use new method of photosystem_coefficients!
#     2022-Feb-28: use updated light_limited_rate! function
#     2022-Feb-28: use updated photosystem_electron_transport! function (twice in thf function)
#     2022-Feb-28: add support to C3CytochromeModel
#     2022-Jun-27: remove apar from input variable list of light_limited_rate!
#     2022-Jun-27: use ClimaCache v0.4, where Leaf.apar is renamed to Leaf.ppar
#     2022-Jun-28: unpack the constant fields
# Bug fixes
#     2022-Jan-24: fix PSM abstraction in colimit_photosynthesis! function
# To do
#     TODO: update leaf T in StomataModels module or higher level
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, g_lc::FT = leaf.g_CO₂) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ diffusive conductance, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, default is `leaf.g_CO₂`

Note that CO₂ partial pressures at leaf surface (stomatal opening) and in leaf internal airspace are updated, and then electron transport is updated again based on this CO₂ partial pressure.

---
# Examples
```julia
leaf = Leaf{Float64}("C3");
air  = AirLayer{Float64}();
mode = GCO₂Mode();
leaf_photosynthesis!(leaf, air, mode);
leaf_photosynthesis!(leaf, air, mode, 0.1);
```
"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, g_lc::FT = leaf.g_CO₂) where {FT<:AbstractFloat} = (
    @unpack PRC, PSM = leaf;

    leaf.g_CO₂ = g_lc;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    # leaf.p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    if leaf.t != leaf._t
        photosystem_temperature_dependence!(PSM, air, leaf.t);
    end;
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf.p_CO₂_i);
    rubisco_limited_rate!(PSM, air, leaf.g_CO₂);
    light_limited_rate!(PSM, PRC, air, leaf.g_CO₂);
    product_limited_rate!(PSM, air, leaf.g_CO₂);
    colimit_photosynthesis!(PSM);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf.p_CO₂_i = air.p_CO₂ - PSM.a_net / leaf.g_CO₂   * air.P_AIR * FT(1e-6);
    leaf.p_CO₂_s = air.p_CO₂ - PSM.a_net / leaf.g_CO₂_b * air.P_AIR * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf.p_CO₂_i);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaf.ppar);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-28: add method for Leaves1D
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::GCO₂Mode) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ diffusive conductance, given
- `leaves` `Leaves1D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates

---
# Examples
```julia
leaves = Leaves1D{Float64}("C3");
air    = AirLayer{Float64}();
mode   = GCO₂Mode();
leaf_photosynthesis!(leaves, air, mode);
```
"""
leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::GCO₂Mode) where {FT<:AbstractFloat} = (
    @unpack PRC, PSM = leaves;

    # leaf.p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    # loop through the leaves.ppar
    for _i in eachindex(leaves.ppar)
        photosystem_temperature_dependence!(PSM, air, leaves.t[_i]);
        photosystem_electron_transport!(PSM, PRC, leaves.ppar[_i], leaves.p_CO₂_i[_i]);
        rubisco_limited_rate!(PSM, air, leaves.g_CO₂[_i]);
        light_limited_rate!(PSM, PRC, air, leaves.g_CO₂[_i]);
        product_limited_rate!(PSM, air, leaves.g_CO₂[_i]);
        colimit_photosynthesis!(PSM);

        # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
        leaves.p_CO₂_i[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂[_i]   * air.P_AIR * FT(1e-6);
        leaves.p_CO₂_s[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b[_i] * air.P_AIR * FT(1e-6);

        # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
        photosystem_electron_transport!(PSM, PRC, leaves.ppar[_i], leaves.p_CO₂_i[_i]);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar[_i]);

        # save the rates and to leaves
        leaves.a_net[_i] = PSM.a_net;
        leaves.a_gross[_i] = PSM.a_gross;
    end;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-28: add method for Leaves2D
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::GCO₂Mode) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ diffusive conductance, given
- `leaves` `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates

---
# Examples
```julia
leaves = Leaves2D{Float64}("C3");
air    = AirLayer{Float64}();
mode   = GCO₂Mode();
leaf_photosynthesis!(leaves, air, mode);
```
"""
leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::GCO₂Mode) where {FT<:AbstractFloat} = (
    @unpack PRC, PSM = leaves;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaves.t != leaves._t
        photosystem_temperature_dependence!(PSM, air, leaves.t);
    end;

    # leaf.p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    # loop through sunlit leaves
    for _i in eachindex(leaves.ppar_sunlit)
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves.p_CO₂_i_sunlit[_i]);
        rubisco_limited_rate!(PSM, air, leaves.g_CO₂_sunlit[_i]);
        light_limited_rate!(PSM, PRC, air, leaves.g_CO₂_sunlit[_i]);
        product_limited_rate!(PSM, air, leaves.g_CO₂_sunlit[_i]);
        colimit_photosynthesis!(PSM);

        # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
        leaves.p_CO₂_i_sunlit[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_sunlit[_i] * air.P_AIR * FT(1e-6);
        leaves.p_CO₂_s_sunlit[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b          * air.P_AIR * FT(1e-6);

        # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves.p_CO₂_i_sunlit[_i]);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar_sunlit[_i]);

        # save the rates and to leaves
        leaves.a_net_sunlit[_i] = PSM.a_net;
        leaves.a_gross_sunlit[_i] = PSM.a_gross;
        leaves.ϕ_f_sunlit[_i] = PRC.ϕ_f;
    end;

    # run the model for shaded leaves
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves.p_CO₂_i_shaded);
    rubisco_limited_rate!(PSM, air, leaves.g_CO₂_shaded);
    light_limited_rate!(PSM, PRC, air, leaves.g_CO₂_shaded);
    product_limited_rate!(PSM, air, leaves.g_CO₂_shaded);
    colimit_photosynthesis!(PSM);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaves.p_CO₂_i_shaded = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_shaded * air.P_AIR * FT(1e-6);
    leaves.p_CO₂_s_shaded = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b      * air.P_AIR * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves.p_CO₂_i_shaded);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaves.ppar_shaded);

    # save the rates and to leaves
    leaves.a_net_shaded = PSM.a_net;
    leaves.a_gross_shaded = PSM.a_gross;
    leaves.ϕ_f_shaded = PRC.ϕ_f;

    return nothing
);


######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for MonoElementSPAC
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `MonoElementSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`
"""
leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat} = (
    @unpack AIR, LEAF = spac;

    leaf_photosynthesis!(LEAF, AIR, mode);

    return nothing
);


######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jun-29: add method for MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, `MonoMLTreeSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`
"""
leaf_photosynthesis!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat} = (
    @unpack AIR, LEAVES, LEAVES_INDEX = spac;

    for _i in eachindex(LEAVES)
        leaf_photosynthesis!(LEAVES[_i], AIR[LEAVES_INDEX[_i]], mode);
    end;

    return nothing
);
