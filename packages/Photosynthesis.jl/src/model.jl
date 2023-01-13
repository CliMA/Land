#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function leaf_photosynthesis!
#
#######################################################################################################################################################################################################
"""
Per refactored Photosynthesis module, the only things one need to know is the public function `leaf_photosynthesis!` and some construtors from `EmeraldNamespace`. See the examples in the methods
    below for details about how to use the function. The steps for computing photosynthetic rates are

- Update temperature dependent variables using [`photosystem_temperature_dependence!`](@ref)
- Calculate electron transport rate using [`photosystem_electron_transport!`](@ref)
- Calculate RubisCO limited rate using [`rubisco_limited_rate!`](@ref)
- Calculate light limited rate using [`light_limited_rate!`](@ref)
- Calculate product limited rate using [`product_limited_rate!`](@ref)
- Calculate gross and net rates using [`colimit_photosynthesis!`](@ref)
- Update fluorescence related variables using [`photosystem_coefficients!`](@ref)

"""
function leaf_photosynthesis! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-07: add method to compute photosynthetic rates only
#     2022-Jul-25: abstractize method to support Leaves1D
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT = lf.t) where {FT<:AbstractFloat}
    leaf_photosynthesis!(lf::Leaves1D{FT}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ partial pressure (for StomataModels.jl temporary use), given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, default is `leaf._g_CO₂`
- `ppar` APAR used for photosynthesis
- `t` Leaf temperature in `[K]`

"""
leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves2D{FT}}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT = lf.t) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = lf;

    photosystem_temperature_dependence!(PSM, air, t);
    photosystem_electron_transport!(PSM, PRC, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(PSM, air, g_lc; β = FT(1));
    light_limited_rate!(PSM, PRC, air, g_lc; β = FT(1));
    product_limited_rate!(PSM, air, g_lc; β = FT(1));
    colimit_photosynthesis!(PSM; β = FT(1));

    return nothing
);

leaf_photosynthesis!(lf::Leaves1D{FT}, air::AirLayer{FT}, g_lc::FT, ppar::FT, t::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = lf;

    photosystem_temperature_dependence!(PSM, air, t);
    photosystem_electron_transport!(PSM, PRC, ppar, FT(20); β = FT(1));
    rubisco_limited_rate!(PSM, air, g_lc; β = FT(1));
    light_limited_rate!(PSM, PRC, air, g_lc; β = FT(1));
    product_limited_rate!(PSM, air, g_lc; β = FT(1));
    colimit_photosynthesis!(PSM; β = FT(1));

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jul-12: add method to account for tuning factor at leaf level
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates based on CO₂ partial pressure or CO₂ conductance, given
- `lf` `Leaf`, `Leaves1D`, or `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` or `PCO₂Mode` that uses CO₂ conductance or partial pressure to compute photosynthetic rates

"""
leaf_photosynthesis!(lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}}, air::AirLayer{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat} = leaf_photosynthesis!(lf, air, mode, lf.SM);

leaf_photosynthesis!(
            lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            sm::AbstractStomataModel{FT}
) where {FT<:AbstractFloat} = leaf_photosynthesis!(lf, air, mode, FT(1));

leaf_photosynthesis!(
            lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            sm::Union{BallBerrySM{FT}, GentineSM{FT}, LeuningSM{FT}, MedlynSM{FT}}
) where {FT<:AbstractFloat} = leaf_photosynthesis!(lf, air, mode, sm.β, sm.β.PARAM_Y);

leaf_photosynthesis!(
            lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            β::BetaFunction{FT},
            param_y::BetaParameterG1
) where {FT<:AbstractFloat} = leaf_photosynthesis!(lf, air, mode, FT(1));

leaf_photosynthesis!(
            lf::Union{Leaf{FT}, Leaves1D{FT}, Leaves2D{FT}},
            air::AirLayer{FT},
            mode::Union{GCO₂Mode, PCO₂Mode},
            β::BetaFunction{FT},
            param_y::BetaParameterVcmax
) where {FT<:AbstractFloat} = leaf_photosynthesis!(lf, air, mode, β.β₁);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: set a default p_i from leaf to combine two methods
#     2022-Jan-14: do not update temperature to avoid its impact on plant hydraulics
#     2022-Jan-14: use colimit function to compute a_gross and a_net
#     2022-Jan-14: set a default g_lc from leaf to combine two methods
#     2022-Jan-18: add p_i to electron transport function input variables
#     2022-Jan-24: fix PSM abstraction in colimit_photosynthesis! function
#     2022-Feb-07: use new method of photosystem_coefficients!
#     2022-Feb-28: use updated light_limited_rate! function
#     2022-Feb-28: add support to C3CytochromeModel
#     2022-Jun-27: remove ppar from input variable list of light_limited_rate!
#     2022-Jun-28: add method for Leaves1D
#     2022-Jun-28: add method for Leaves2D
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-07: save a_net and a_gross to Leaf (as PSM may be used for temporary calculations)
#     2022-Jul-12: use β as a must have option (and thus this function becomes a core function of the one above)
#     2022-Jul-28: move temperature control to function photosystem_temperature_dependence!
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT) where {FT<:AbstractFloat}
    leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT) where {FT<:AbstractFloat}
    leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT) where {FT<:AbstractFloat}
    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT) where {FT<:AbstractFloat}
    leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT) where {FT<:AbstractFloat}
    leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates (this method not meant for public usage, use it with caution), given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `leaves` `Leaves1D` or `Leaves2D` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `GCO₂Mode` or `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = leaf;

    photosystem_temperature_dependence!(PSM, air, leaf.t);
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf._p_CO₂_i; β = β);
    rubisco_limited_rate!(PSM, leaf._p_CO₂_i; β = β);
    light_limited_rate!(PSM);
    product_limited_rate!(PSM, leaf._p_CO₂_i; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaf.ppar; β = β);

    # save the rates and to leaf
    leaf.a_net = PSM.a_net;
    leaf.a_gross = PSM.a_gross;

    return nothing
);

leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = leaves;

    # loop through the ppars
    for _i in eachindex(leaves.ppar)
        # update TD parameters everytime for sunlit and shaded leaves
        photosystem_temperature_dependence!(PSM, air, leaves.t[_i]);

        # calculate the photosynthetic rates
        photosystem_electron_transport!(PSM, PRC, leaves.ppar[_i], leaves._p_CO₂_i[_i]; β = β);
        rubisco_limited_rate!(PSM, leaves._p_CO₂_i[_i]; β = β);
        light_limited_rate!(PSM);
        product_limited_rate!(PSM, leaves._p_CO₂_i[_i]; β = β);
        colimit_photosynthesis!(PSM; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar[_i]; β = β);

        # save the rates and to leaves
        leaves.a_net[_i] = PSM.a_net;
        leaves.a_gross[_i] = PSM.a_gross;
    end;

    return nothing
);

leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::PCO₂Mode, β::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = leaves;

    photosystem_temperature_dependence!(PSM, air, leaves.t);

    # loop through the ppars for sunlit leaves
    for _i in eachindex(leaves.ppar_sunlit)
        # calculate the photosynthetic rates
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves._p_CO₂_i_sunlit[_i]; β = β);
        rubisco_limited_rate!(PSM, leaves._p_CO₂_i_sunlit[_i]; β = β);
        light_limited_rate!(PSM);
        product_limited_rate!(PSM, leaves._p_CO₂_i_sunlit[_i]; β = β);
        colimit_photosynthesis!(PSM; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar_sunlit[_i]; β = β);

        # save the rates and to leaves
        leaves.a_net_sunlit[_i] = PSM.a_net;
        leaves.a_gross_sunlit[_i] = PSM.a_gross;
        leaves.ϕ_f_sunlit[_i] = PRC.ϕ_f;
    end;

    # run the model for shaded leaves
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves._p_CO₂_i_shaded; β = β);
    rubisco_limited_rate!(PSM, leaves._p_CO₂_i_shaded; β = β);
    light_limited_rate!(PSM);
    product_limited_rate!(PSM, leaves._p_CO₂_i_shaded; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaves.ppar_shaded; β = β);

    # save the rates and to leaves
    leaves.a_net_shaded = PSM.a_net;
    leaves.a_gross_shaded = PSM.a_gross;
    leaves.ϕ_f_shaded = PRC.ϕ_f;

    return nothing
);

leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = leaf;

    # leaf._p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    photosystem_temperature_dependence!(PSM, air, leaf.t);
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf._p_CO₂_i; β = β);
    rubisco_limited_rate!(PSM, air, leaf._g_CO₂; β = β);
    light_limited_rate!(PSM, PRC, air, leaf._g_CO₂; β = β);
    product_limited_rate!(PSM, air, leaf._g_CO₂; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaf._p_CO₂_i = air.p_CO₂ - PSM.a_net / leaf._g_CO₂   * air.P_AIR * FT(1e-6);
    leaf._p_CO₂_s = air.p_CO₂ - PSM.a_net / leaf.g_CO₂_b * air.P_AIR * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(PSM, PRC, leaf.ppar, leaf._p_CO₂_i; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaf.ppar; β = β);

    # save the rates and to leaf
    leaf.a_net = PSM.a_net;
    leaf.a_gross = PSM.a_gross;

    return nothing
);

leaf_photosynthesis!(leaves::Leaves1D{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = leaves;

    # leaf._p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    # loop through the leaves.ppar
    for _i in eachindex(leaves.ppar)
        photosystem_temperature_dependence!(PSM, air, leaves.t[_i]);
        photosystem_electron_transport!(PSM, PRC, leaves.ppar[_i], leaves._p_CO₂_i[_i]; β = β);
        rubisco_limited_rate!(PSM, air, leaves._g_CO₂[_i]; β = β);
        light_limited_rate!(PSM, PRC, air, leaves._g_CO₂[_i]; β = β);
        product_limited_rate!(PSM, air, leaves._g_CO₂[_i]; β = β);
        colimit_photosynthesis!(PSM; β = β);

        # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
        leaves._p_CO₂_i[_i] = air.p_CO₂ - PSM.a_net / leaves._g_CO₂[_i]   * air.P_AIR * FT(1e-6);
        leaves._p_CO₂_s[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b[_i] * air.P_AIR * FT(1e-6);

        # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
        photosystem_electron_transport!(PSM, PRC, leaves.ppar[_i], leaves._p_CO₂_i[_i]; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar[_i]; β = β);

        # save the rates and to leaves
        leaves.a_net[_i] = PSM.a_net;
        leaves.a_gross[_i] = PSM.a_gross;
    end;

    return nothing
);

leaf_photosynthesis!(leaves::Leaves2D{FT}, air::AirLayer{FT}, mode::GCO₂Mode, β::FT) where {FT<:AbstractFloat} = (
    (; PRC, PSM) = leaves;

    photosystem_temperature_dependence!(PSM, air, leaves.t);

    # leaf._p_CO₂_i is not accurate here in the first call, thus need a second call after p_CO₂_i is analytically resolved
    # loop through sunlit leaves
    for _i in eachindex(leaves.ppar_sunlit)
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves._p_CO₂_i_sunlit[_i]; β = β);
        rubisco_limited_rate!(PSM, air, leaves._g_CO₂_sunlit[_i]; β = β);
        light_limited_rate!(PSM, PRC, air, leaves._g_CO₂_sunlit[_i]; β = β);
        product_limited_rate!(PSM, air, leaves._g_CO₂_sunlit[_i]; β = β);
        colimit_photosynthesis!(PSM; β = β);

        # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
        leaves._p_CO₂_i_sunlit[_i] = air.p_CO₂ - PSM.a_net / leaves._g_CO₂_sunlit[_i] * air.P_AIR * FT(1e-6);
        leaves._p_CO₂_s_sunlit[_i] = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b          * air.P_AIR * FT(1e-6);

        # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
        photosystem_electron_transport!(PSM, PRC, leaves.ppar_sunlit[_i], leaves._p_CO₂_i_sunlit[_i]; β = β);

        # update the fluorescence related parameters
        photosystem_coefficients!(PSM, PRC, leaves.ppar_sunlit[_i]; β = β);

        # save the rates and to leaves
        leaves.a_net_sunlit[_i] = PSM.a_net;
        leaves.a_gross_sunlit[_i] = PSM.a_gross;
        leaves.ϕ_f_sunlit[_i] = PRC.ϕ_f;
    end;

    # run the model for shaded leaves
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves._p_CO₂_i_shaded; β = β);
    rubisco_limited_rate!(PSM, air, leaves._g_CO₂_shaded; β = β);
    light_limited_rate!(PSM, PRC, air, leaves._g_CO₂_shaded; β = β);
    product_limited_rate!(PSM, air, leaves._g_CO₂_shaded; β = β);
    colimit_photosynthesis!(PSM; β = β);

    # update CO₂ partial pressures at the leaf surface and internal airspace (evaporative front)
    leaves._p_CO₂_i_shaded = air.p_CO₂ - PSM.a_net / leaves._g_CO₂_shaded * air.P_AIR * FT(1e-6);
    leaves._p_CO₂_s_shaded = air.p_CO₂ - PSM.a_net / leaves.g_CO₂_b      * air.P_AIR * FT(1e-6);

    # update leaf ETR again to ensure that j_pot and e_to_c are correct for C3CytochromeModel
    photosystem_electron_transport!(PSM, PRC, leaves.ppar_shaded, leaves._p_CO₂_i_shaded; β = β);

    # update the fluorescence related parameters
    photosystem_coefficients!(PSM, PRC, leaves.ppar_shaded; β = β);

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
#     2022-Jun-29: add method for MonoMLGrassSPAC, MonoMLPalmSPAC, MonoMLTreeSPAC
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-13: redirect the wrapper function to the method at leaf level
#
#######################################################################################################################################################################################################
"""

    leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}
    leaf_photosynthesis!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates for SPAC, given
- `spac` `MonoElementSPAC`, `MonoMLGrassSPAC`, `MonoMLPalmSPAC`, or `MonoMLTreeSPAC` type SPAC
- `mode` `GCO₂Mode` or `PCO₂Mode`

"""
leaf_photosynthesis!(spac::MonoElementSPAC{FT}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat} = leaf_photosynthesis!(spac.LEAF, spac.AIR, mode);

leaf_photosynthesis!(spac::Union{MonoMLGrassSPAC{FT}, MonoMLPalmSPAC{FT}, MonoMLTreeSPAC{FT}}, mode::Union{GCO₂Mode, PCO₂Mode}) where {FT<:AbstractFloat} = (
    (; AIR, LEAVES, LEAVES_INDEX) = spac;

    for _i in eachindex(LEAVES)
        leaf_photosynthesis!(LEAVES[_i], AIR[LEAVES_INDEX[_i]], mode);
    end;

    return nothing
);
