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
    _p = air.p_CO₂;
    _r = ps.r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _b - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    ps.a_c = _an + _r;

    return nothing
);


"""
    leaf_ac!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the RubisCO limited photosynthetic rate in conductance mode, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in the method
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this methid
"""
leaf_ac!(ps::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
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
    _p = air.p_CO₂;
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
    leaf_aj!(ps::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `ps` `C4VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in this methid
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this methid
"""
leaf_aj!(ps::C4VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
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
    leaf_ap!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat}

Update the electron transport limited photosynthetic rate in conductance mode, given
- `ps` `C3VJPModel` structure for C3 photosynthesis model
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure, not used in this method
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, not used in this method
"""
leaf_ap!(ps::C3VJPModel{FT}, air::AirLayer{FT}, g_lc::FT) where {FT<:AbstractFloat} = (
    ps.a_p = ps.v_cmax / 2;

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
    _p = air.p_CO₂;
    _r = ps.r_d;

    _qa = _f;
    _qb = _f*_r - _p - _d - _a*_f;
    _qc = _a*_p - _r*(_p + _d);
    _an = lower_quadratic(_qa, _qb, _qc);

    ps.a_p = _an + _r;

    return nothing
);


"""
This function calculates leaf photosynthetic rates and save the values within the structure. Supported methods are

$(METHODLIST)
"""
function leaf_photosynthesis! end


"""
    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `p_i` Internal CO₂ partial pressure in `Pa`, default is `leaf.p_CO₂_i`
"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat} = (
    leaf.p_CO₂_i = p_i;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaf.t != leaf._t
        photosystem_temperature_dependence!(leaf.PSM, air, leaf.t);
    end;
    photosystem_electron_transport!(leaf.PSM, leaf.PRC, leaf.apar);
    leaf_ac!(leaf.PSM, leaf.p_CO₂_i);
    leaf_aj!(leaf.PSM, leaf.p_CO₂_i);
    leaf_ap!(leaf.PSM, leaf.p_CO₂_i);

    # TODO: add a function leaf_ag! to use colimit functions
    leaf.PSM.a_gross = min(leaf.PSM.a_c, leaf.PSM.a_j, leaf.PSM.a_p);
    leaf.PSM.a_net   = leaf.PSM.a_gross - leaf.PSM.r_d;

    return nothing
);


"""
    leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::PCO₂Mode, p_i::FT = leaf.p_CO₂_i) where {FT<:AbstractFloat}

Updates leaf photosynthetic rates, given
- `leaf` `Leaf` type structure that stores biophysical, reaction center, and photosynthesis model structures
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `mode` `PCO₂Mode` that uses CO₂ partial pressure to compute photosynthetic rates
- `g_lc` Leaf diffusive conductance to CO₂ in `[mol m⁻² s⁻¹]`, default is `leaf.g_CO₂`
"""
leaf_photosynthesis!(leaf::Leaf{FT}, air::AirLayer{FT}, mode::GCO₂Mode, g_lc::FT = leaf.g_CO₂) where {FT<:AbstractFloat} = (
    leaf.g_CO₂ = g_lc;

    # because xylem parameters and vapor pressure are also temperature dependent, do not change leaf._t here!
    if leaf.t != leaf._t
        photosystem_temperature_dependence!(leaf.PSM, air, leaf.t);
    end;
    photosystem_electron_transport!(leaf.PSM, leaf.PRC, leaf.apar);
    leaf_ac!(leaf.PSM, air, leaf.g_CO₂);
    leaf_aj!(leaf.PSM, air, leaf.g_CO₂);
    leaf_ap!(leaf.PSM, air, leaf.g_CO₂);

    # TODO: add a function leaf_ag! to use colimit functions
    leaf.PSM.a_gross = min(leaf.PSM.a_c, leaf.PSM.a_j, leaf.PSM.a_p);
    leaf.PSM.a_net   = leaf.PSM.a_gross - leaf.PSM.r_d;

    leaf.p_CO₂_i = air.p_CO₂ - leaf.PSM.a_net / leaf.g_CO₂   * air.P_AIR * FT(1e-6);
    leaf.p_CO₂_s = air.p_CO₂ - leaf.PSM.a_net / leaf.g_CO₂_b * air.P_AIR * FT(1e-6);

    return nothing
);
