#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function product_limited_rate!
#
#######################################################################################################################################################################################################
"""
This function updates the product limited photosynthetic rate. Supported methods are

$(METHODLIST)

"""
function product_limited_rate! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add input variable p_i to make the code more modular
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the product limited photosynthetic rate, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
product_limited_rate!(psm::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    psm.a_p = psm.v_cmax / 2;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-18: add support to C3CytochromeModel
#     2022-Jan-18: add input variable p_i to make the code more modular
#     2022-Jan-24: fix documentation
# Bug fix
#     2022-Jan-18: j_psi_p = j_psii_p * η (was j_psii_c * η)
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::C3CytochromeModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the product limited photosynthetic rate, given
- `psm` `C3CytochromeModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
product_limited_rate!(psm::C3CytochromeModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;

    psm.a_p      = psm.v_cmax / 2;
    psm.j_psii_p = psm.a_c * (EFF_1*p_i + EFF_2*psm.γ_star) / (p_i - psm.γ_star);
    psm.j_psi_p  = psm.j_psii_p* psm.η;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add input variable p_i to make the code more modular
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the product limited photosynthetic rate, given
- `psm` `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
product_limited_rate!(psm::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    psm.a_p = psm.v_pmax * p_i / (p_i + psm.k_pep);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add this new method to simplify the multiple dispatch of leaf_photosynthesis!
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in this method
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this method
"""
product_limited_rate!(psm::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    psm.a_p = psm.v_cmax / 2;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add input variable g_lc to make the code more modular
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C4VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
product_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    _a = psm.v_pmax;
    _d = psm.k_pep;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = psm.r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    psm.a_p = _an + _r;

    return nothing
);
