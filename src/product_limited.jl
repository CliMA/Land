#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function product_limited_rate!
#     2022-Jul-13: deflate documentation
#
#######################################################################################################################################################################################################
"""
This function supports two types of calculations:
- Calculate the rate from internal CO₂
- Calculate the rate from CO₂ conductance by solving a quadratic function

"""
function product_limited_rate! end


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Jan-14: add input variable p_i to make the code more modular
#     2022-Jan-18: add support to C3CytochromeModel
#     2022-Jan-24: fix documentation
#     2022-Feb-07: add C3CytochromeModel support into Union
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-13: deflate documentation
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    product_limited_rate!(psm::C4VJPModel{FT}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Update the product limited photosynthetic rate, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (psm._a_p = β * psm._v_cmax / 2; return nothing);

product_limited_rate!(psm::C4VJPModel{FT}, p_i::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (psm._a_p = β * psm._v_pmax * p_i / (p_i + psm._k_pep); return nothing);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add this new method to simplify the multiple dispatch of leaf_photosynthesis!
#     2022-Jan-14: add input variable g_lc to make the code more modular
#     2022-Jan-24: fix documentation
#     2022-Feb-28: add C3CytochromeModel support
#     2022-Jul-01: add β to variable list to account for Vmax downregulation used in CLM5
#     2022-Jul-13: deflate documentation
#
#######################################################################################################################################################################################################
"""

    product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}
    product_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C3CytochromeModel`, `C3VJPModel`, or `C4VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in this method
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this method
- `β` Tuning factor to downregulate effective Vmax, Jmax, and Rd

"""
product_limited_rate!(psm::Union{C3CytochromeModel{FT}, C3VJPModel{FT}}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (psm._a_p = β * psm._v_cmax / 2; return nothing);

product_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT; β::FT = FT(1)) where {FT<:AbstractFloat} = (
    _a = β * psm._v_pmax;
    _d = psm._k_pep;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = β * psm._r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    psm._a_p = _an + _r;

    return nothing
);
