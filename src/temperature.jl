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


"""
    temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the correction ratio for a temperature dependent variable, given
- `td` `Arrhenius` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

---
# Examples
```julia
tdep = Arrhenius{Float64}(298.15, 28208.88, 36380.0);
cor1 = temperature_correction(tdep, 300.0);
cor2 = temperature_correction(tdep, 300.0; t_ref=290.0);
```
"""
temperature_correction(td::Arrhenius{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat} = exp( td.ΔHA / GAS_R(FT) * (1/t_ref - 1/t) );


"""
    temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the correction ratio for a temperature dependent variable, given
- `td` `ArrheniusPeak` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

---
# Examples
```julia
tdep = Arrhenius{Float64}(298.15, 1, 94800.0, 73300.0, 250.0);
cor1 = temperature_correction(tdep, 300.0);
cor2 = temperature_correction(tdep, 300.0; t_ref=290.0);
```
"""
temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref = td.T_REF) where {FT<:AbstractFloat} = (
    @unpack ΔHA, ΔHD, ΔSV = td;

    # _f_a: activation correction, _f_b: de-activation correction
    _f_a = exp( ΔHA / GAS_R(FT) * (1 / t_ref - 1 / t) );
    _f_b = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t_ref))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * t)));

    return _f_a * _f_b
);


"""
    temperature_correction(td::ArrheniusPeak{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the correction ratio for a temperature dependent variable, given
- `td` `Q10` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

---
# Examples
```julia
tdep = Q10{Float64}(298.15, 1.0, 1.3);
cor1 = temperature_correction(tdep, 300.0);
cor2 = temperature_correction(tdep, 300.0; t_ref=290.0);
```
"""
temperature_correction(td::Q10{FT}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat} = td.Q_10 ^ ( (t - t_ref) / 10 );


"""
    temperature_corrected_value(td::Union{Arrhenius{FT}, ArrheniusPeak{FT}, Q10{FT}}, t::FT; t_ref::FT = td.T_REF) where {FT<:AbstractFloat}

Return the temperature corrected value, given
- `td` `Q10` type temperature dependency struture
- `t` Target temperature in `K`
- `t_ref` Reference temperature in `K`, default is `td.T_REF` (298.15 K)

---
# Examples
```julia
tdep = Q10{Float64}(298.15, 1.0, 1.3);
val1 = temperature_corrected_value(tdep, 300.0);
val2 = temperature_corrected_value(tdep, 300.0; t_ref=290.0);
```
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


"""
    photosystem_temperature_dependence!(ps::C3VJPSystem{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}

Update the temperature dependencies of C3 photosynthesis system, given
- `ps` `C3VJPSystem` structure for C3 photosynthesis system
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`

---
# Examples
```julia
ps = C3VJPSystem{Float64}();
air = AirLayer{Float64}();
photosystem_temperature_dependence!(ps, air, 300.0);
```
"""
photosystem_temperature_dependence!(ps::C3VJPSystem{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat} = (
    ps.r_d    = ps.r_d25    * temperature_correction(ps.TD_R, t);
    ps.v_cmax = ps.v_cmax25 * temperature_correction(ps.TD_VCMAX, t);
    ps.j_max  = ps.j_max25  * temperature_correction(ps.TD_JMAX, t);
    ps.k_c    = temperature_corrected_value(ps.TD_KC, t);
    ps.k_o    = temperature_corrected_value(ps.TD_KO, t);
    ps.γ_star = temperature_corrected_value(ps.TD_Γ, t);
    ps.k_m    = ps.k_c * (1 + air.P_O2 / ps.k_o);

    return nothing
);


"""
    photosystem_temperature_dependence!(ps::C4VJPSystem{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat}

Update the temperature dependencies of C3 photosynthesis system, given
- `ps` `C4VJPSystem` structure for C3 photosynthesis system
- `air` `AirLayer` structure for environmental conditions like O₂ partial pressure
- `t` Target temperature in `K`

---
# Examples
```julia
ps = C4VJPSystem{Float64}();
air = AirLayer{Float64}();
photosystem_temperature_dependence!(ps, air, 300.0);
```
"""
photosystem_temperature_dependence!(ps::C4VJPSystem{FT}, air::AirLayer{FT}, t::FT) where {FT<:AbstractFloat} = (
    ps.r_d    = ps.r_d25    * temperature_correction(ps.TD_R, t);
    ps.v_cmax = ps.v_cmax25 * temperature_correction(ps.TD_VCMAX, t);
    ps.v_pmax = ps.v_pmax25 * temperature_correction(ps.TD_VPMAX, t);
    ps.k_pep  = temperature_corrected_value(ps.TD_KPEP, t);

    return nothing
);
