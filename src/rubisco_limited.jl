#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: refactor the function rubisco_limited_rate!
#
#######################################################################################################################################################################################################
"""
This function updates the RubisCO limited photosynthetic rate. Supported methods are

$(METHODLIST)

"""
function rubisco_limited_rate! end


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add input variable p_i to make the code more modular
#
#######################################################################################################################################################################################################
"""
    rubisco_limited_rate!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
rubisco_limited_rate!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.a_c = ps.v_cmax * (p_i - ps.γ_star) / (p_i + ps.k_m);

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add input variable p_i to make the code more modular
#
#######################################################################################################################################################################################################
"""
    rubisco_limited_rate!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate, given
- `ps` `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
rubisco_limited_rate!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.a_c = ps.v_cmax;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add input variable g_lc to make the code more modular
#
#######################################################################################################################################################################################################
"""
    rubisco_limited_rate!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate in conductance mode, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
rubisco_limited_rate!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    _a = ps.v_max;
    _b = ps.v_max * ps.γ_star;
    _d = ps.k_m;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_CO₂;
    _r = ps.r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    ps.a_c = _an + _r;

    return nothing
);


#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Jan-14: add this new method to simplify the multiple dispatch of leaf_photosynthesis!
#
#######################################################################################################################################################################################################
"""
    rubisco_limited_rate!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate in conductance mode, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in the method
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this methid
"""
rubisco_limited_rate!(ps::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    ps.a_c = ps.v_cmax;

    return nothing
);
