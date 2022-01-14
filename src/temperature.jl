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
