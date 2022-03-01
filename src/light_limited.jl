#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function light_limited_rate!
#
#######################################################################################################################################################################################################
"""
This function updates the electron transport limited photosynthetic rate. Supported methods are

$(METHODLIST)

"""
function light_limited_rate! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: unpack CONSTANT from the input variables only
#     2022-Jan-14: add p_i to input list to make the code more modular
#     2022-Jan-24: add C3CytochromeModel support in a Union
#     2022-Jan-24: fix documentation
#     2022-Feb-07: move C3CytochromeModel support out given different field name
#     2022-Feb-28: move e_to_c calculation to photosystem_electron_transport!
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3VJPModel{FT}) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
"""
light_limited_rate!(psm::C3VJPModel{FT}) where {FT<:AbstractFloat} = (
    psm.a_j = psm.j * psm.e_to_c;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: add p_i to input list to make the code more modular
#     2022-Jan-24: fix documentation
#     2022-Feb-07: remove duplicated j (using j_pot is enough)
#     2022-Feb-07: move C3CytochromeModel support out given different field name
#     2022-Feb-28: move e_to_c calculation to photosystem_electron_transport!
#     2022-Feb-28: add C4VJPModel as a Union type of psm
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `psm` `C3CytochromeModel` or `C4VJPModel` structure for C3/C4 photosynthesis model
"""
light_limited_rate!(psm::Union{C3CytochromeModel{FT}, C4VJPModel{FT}}) where {FT<:AbstractFloat} = (
    psm.a_j = psm.j_pot * psm.e_to_c;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Feb-28: add C3CytochromeModel support
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C3CytochromeModel` structure for C3 photosynthesis model
- `rc` `CytochromeReactionCenter` type photosynthesis system reaction center
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
light_limited_rate!(psm::C3CytochromeModel{FT}, rc::CytochromeReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;
    @unpack F_PSI, Η_C, Η_L, Φ_PSI_MAX = rc;

    _j_psi = psm.v_qmax * apar * F_PSI * Φ_PSI_MAX / (psm.v_qmax + apar * F_PSI * Φ_PSI_MAX);
    _eff_a = 1 - Η_L / Η_C;
    _eff_b = 1 / Η_C;
    _eff_1 = _eff_a * EFF_1 + 3 * _eff_b;
    _eff_2 = _eff_a * EFF_2 + 7 * _eff_b;

    _a = _j_psi;
    _b = _j_psi * psm.γ_star;
    _c = _eff_1;
    _d = _eff_2 * psm.γ_star;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = psm.r_d;

    _qa = _c * _f;
    _qb = _c*_f*_r - _c*_p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_c*_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    psm.a_j = _an + _r;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add g_lc to input list to make the code more modular
#     2022-Jan-24: fix a bug in field name e_to_c in psm
#     2022-Jan-24: fix documentation
#     2022-Feb-28: move e_to_c calculation to photosystem_electron_transport!
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center, not used in this methid
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`, not used in this methid
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
light_limited_rate!(psm::C3VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat} = (
    _a = psm.j;
    _b = psm.j * psm.γ_star;
    _c = psm.EFF_1;
    _d = psm.EFF_2 * psm.γ_star;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = psm.r_d;

    _qa = _c * _f;
    _qb = _c*_f*_r - _c*_p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_c*_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    psm.a_j = _an + _r;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add this new method to simplify the multiple dispatch of leaf_photosynthesis!
#     2022-Jan-24: fix documentation
#     2022-Feb-07: remove duplicated j (using j_pot is enough)
#     2022-Feb-28: move e_to_c calculation to photosystem_electron_transport!
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C4VJPModel` structure for C3 photosynthesis model
- `rc` `VJPReactionCenter` type photosynthesis system reaction center, not used in this methid
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in this methid
- `apar` Absorbed photosynthetically active radiation in `μmol m⁻² s⁻¹`, not used in this methid
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this methid
"""
light_limited_rate!(psm::C4VJPModel{FT}, rc::VJPReactionCenter{FT}, air::AirLayer{FT}, apar::FT, g_lc::FT) where {FT<:AbstractFloat} = (
    psm.a_j = psm.j_pot * psm.e_to_c;

    return nothing
);
