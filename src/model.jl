"""
This function updates the RubisCO limited photosynthetic rate. Supported methods are

$(METHODLIST)
"""
function leaf_ac! end


"""
    leaf_ac!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
leaf_ac!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.a_c = ps.v_cmax * (p_i - ps.γ_star) / (p_i + ps.k_m);

    return nothing
);


"""
    leaf_ac!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate in conductance mode, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
leaf_ac!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    _a = ps.v_max;
    _b = ps.v_max * ps.γ_star;
    _d = ps.k_m;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_co2;
    _r = ps.r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    ps.a_c = _an + _r;

    return nothing
);


"""
    leaf_ac!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate, given
- `ps` `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
leaf_ac!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.a_c = ps.v_cmax;

    return nothing
);


"""
This function updates the electron transport limited photosynthetic rate. Supported methods are

$(METHODLIST)
"""
function leaf_aj! end


"""
    leaf_aj!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
leaf_aj!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.e_to_c = (p_i - ps.γ_star) / (ps.EFF_1 * p_i + ps.EFF_2 * ps.γ_star);
    ps.a_j = ps.j * ps.e_to_c;

    return nothing
);


"""
    leaf_aj!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
leaf_aj!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    _a = ps.j;
    _b = ps.j * ps.γ_star;
    _c = ps.EFF_1;
    _d = ps.EFF_2 * ps.γ_star;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_co2;
    _r = ps.r_d;

    _qa = _c * _f;
    _qb = _c*_f*_r - _c*_p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_c*_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    ps.a_j = _an + _r;
    ps.e2c = ps.a_j / ps.j;

    return nothing
);


"""
    leaf_aj!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate, given
- `ps` `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
leaf_aj!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.e_to_c = 1 / 6;
    ps.a_j = ps.j * ps.e_to_c;

    return nothing
);


"""
This function updates the product limited photosynthetic rate. Supported methods are

$(METHODLIST)
"""
function leaf_ap! end


"""
    leaf_aj!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the product limited photosynthetic rate, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`, not used in this method
"""
leaf_ap!(ps::C3VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.a_p = ps.v_cmax / 2;

    return nothing
);


"""
    leaf_aj!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat}

Update the product limited photosynthetic rate, given
- `ps` `C4VJPModel` structure for C3 photosynthesis model
- `p_i` Internal CO₂ partial pressure in `Pa`
"""
leaf_ap!(ps::C4VJPModel{FT}, p_i::FT) where {FT<:AbstractFloat} = (
    ps.a_p = ps.v_pmax * p_i / (p_i + ps.k_pep);

    return nothing
);


"""
    leaf_ap!(ps::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `ps` `C4VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`
"""
leaf_ap!(ps::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    _a = ps.v_pmax;
    _d = ps.k_pep;
    _f = air.P_AIR / g_lc * FT(1e-6);
    _p = air.p_co2;
    _r = ps.r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    ps.a_p = _an + _r;

    return nothing
);
