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
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
light_limited_rate!(psm::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;

    psm.e_to_c = (p_i - psm.γ_star) / (EFF_1*p_i + EFF_2*psm.γ_star);
    psm.a_j    = psm.j * psm.e_to_c;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this method
# General
#     2022-Feb-07: move C3CytochromeModel support out given different field name
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3CytochromeModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `psm` `C3CytochromeModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
light_limited_rate!(psm::C3CytochromeModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    @unpack EFF_1, EFF_2 = psm;

    psm.e_to_c = (p_i - psm.γ_star) / (EFF_1*p_i + EFF_2*psm.γ_star);
    psm.a_j    = psm.j_pot * psm.e_to_c;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add p_i to input list to make the code more modular
#     2022-Jan-24: fix documentation
#     2022-Feb-07: remove duplicated j (using j_pot is enough)
# To-do
#     TODO: move 1/6 to C4VJPModel definition
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `psm` `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
light_limited_rate!(psm::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    psm.e_to_c = 1 / 6;
    psm.a_j    = psm.j_pot * psm.e_to_c;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add g_lc to input list to make the code more modular
#     2022-Jan-24: fix a bug in field name e_to_c in psm
#     2022-Jan-24: fix documentation
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
light_limited_rate!(psm::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
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

    psm.a_j    = _an + _r;
    psm.e_to_c = psm.a_j / psm.j;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add this new method to simplify the multiple dispatch of leaf_photosynthesis!
#     2022-Jan-24: fix documentation
#     2022-Feb-07: remove duplicated j (using j_pot is enough)
# To-do
#     TODO: move 1/6 to C4VJPModel definition
#
#######################################################################################################################################################################################################
"""

    light_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `psm` `C4VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in this methid
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this methid
"""
light_limited_rate!(psm::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    psm.e_to_c = 1 / 6;
    psm.a_j    = psm.j_pot * psm.e_to_c;

    return nothing
);
