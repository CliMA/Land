#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#
#######################################################################################################################################################################################################
"""
This function calculates the temperature correction ratio relative to the reference temperature. Supported methods are

$(METHODLIST)

"""
function temperature_correction end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-13: add optional input t_ref to allow for manually setting reference temperature
#     2022-Jan-14: remove examples from doc as this function is not meant to be public
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the correction ratio for a temperature dependent variable, given
- `td` `Arrhenius` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)
"""
temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat} = exp( td.ΔHA / GAS_R(FT) * (1/t_ref - 1/t) );


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-13: add optional input t_ref to allow for manually setting reference temperature
#     2022-Jan-14: remove examples from doc as this function is not meant to be public
#     2022-Jan-24: add FT control to r_ref
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the correction ratio for a temperature dependent variable, given
- `td` `ArrheniusPeak` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)
"""
temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat} = (
    @unpack ΔHA, ΔHD, ΔSV = td;

    # _f_a: activation correction, _f_b: de-activation correction
    _f_a = exp( ΔHA / GAS_R(FT) * (1 / t_ref - 1 / t) );
    _f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return _f_a * _f_b
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-13: add optional input t_ref to allow for manually setting reference temperature
#     2022-Jan-14: remove examples from doc as this function is not meant to be public
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the correction ratio for a temperature dependent variable, given
- `td` `Q10` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)
"""
temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat} = td.Q_10 ^ ( (t - t_ref) / 10 );


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-13: add optional input t_ref to allow for manually setting reference temperature
#     2022-Jan-14: remove examples from doc as this function is not meant to be public
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    temperature_corrected_value(td::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the temperature corrected value, given
- `td` `Q10` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)
"""
function temperature_corrected_value(td::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}
    return td.VAL_REF * temperature_correction(td, t; t_ref=t_ref)
end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: use ClimaCache types, which saves photosystem temperature dependencies within the struct
#
#######################################################################################################################################################################################################
"""
This function update all the temperature dependencies of variables in a photosystem, such as Vcmax and Jmax. Supported methods are

$(METHODLIST)

"""
function photosystem_temperature_dependence! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-14: remove examples from doc as this function is not meant to be public
#     2022-Jan-24: fix documentation
#     2022-Mar-03: rearrange the order to make it look better :)
#
#######################################################################################################################################################################################################
"""

    photosystem_temperature_dependence!(psm::C3VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}

Update the temperature dependencies of C3 photosynthesis model, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`
"""
photosystem_temperature_dependence!(psm::C3VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat} = (
    psm.k_c    = temperature_corrected_value(psm.TD_KC, t);
    psm.k_o    = temperature_corrected_value(psm.TD_KO, t);
    psm.γ_star = temperature_corrected_value(psm.TD_Γ , t);
    psm.r_d    = psm.r_d25    * temperature_correction(psm.TD_R, t);
    psm.v_cmax = psm.v_cmax25 * temperature_correction(psm.TD_VCMAX, t);
    psm.j_max  = psm.j_max25  * temperature_correction(psm.TD_JMAX, t);
    psm.k_m    = psm.k_c * (1 + air.P_O₂ / psm.k_o);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Feb-07: add method for C3CytochromeModel photosynthesis model
#     2022-Feb-07: add v_qmax without temperature dependency
#     2022-Mar-01: add temperature dependencies for k_q, v_qmax, η_c, and η_l
#
#######################################################################################################################################################################################################
"""

    photosystem_temperature_dependence!(psm::C3CytochromeModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}

Update the temperature dependencies of C3 photosynthesis model, given
- `psm` `C3CytochromeModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`
"""
photosystem_temperature_dependence!(psm::C3CytochromeModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat} = (
    psm.k_c    = temperature_corrected_value(psm.TD_KC, t);
    psm.k_o    = temperature_corrected_value(psm.TD_KO, t);
    psm.k_q    = temperature_corrected_value(psm.TD_KQ, t);
    psm.γ_star = temperature_corrected_value(psm.TD_Γ , t);
    psm.η_c    = temperature_corrected_value(psm.TD_ΗC, t);
    psm.η_l    = temperature_corrected_value(psm.TD_ΗL, t);
    psm.r_d    = psm.r_d25    * temperature_correction(psm.TD_R, t);
    psm.v_cmax = psm.v_cmax25 * temperature_correction(psm.TD_VCMAX, t);
    psm.k_m    = psm.k_c * (1 + air.P_O₂ / psm.k_o);
    psm.v_qmax = psm.b₆f * psm.k_q;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-13: use ClimaCache types, which uses ΔHA, ΔHD, and ΔSV directly
#     2022-Jan-14: remove examples from doc as this function is not meant to be public
#     2022-Jan-24: fix documentation
#     2022-Mar-03: rearrange the order to make it look better :)
#
#######################################################################################################################################################################################################
"""

    photosystem_temperature_dependence!(psm::C4VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}

Update the temperature dependencies of C3 photosynthesis model, given
- `psm` `C4VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`
"""
photosystem_temperature_dependence!(psm::C4VJPModel{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat} = (
    psm.k_pep  = temperature_corrected_value(psm.TD_KPEP, t);
    psm.r_d    = psm.r_d25    * temperature_correction(psm.TD_R, t);
    psm.v_cmax = psm.v_cmax25 * temperature_correction(psm.TD_VCMAX, t);
    psm.v_pmax = psm.v_pmax25 * temperature_correction(psm.TD_VPMAX, t);

    return nothing
);
