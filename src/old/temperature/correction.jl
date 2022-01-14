###############################################################################
#
# Arrhenius corrections
#
###############################################################################
"""
    temperature_correction(
                td_set::AbstractTDParameterSet{FT},
                T::FT
    ) where {FT<:AbstractFloat}

A correction factor based on arrhenius's fitting procedure, given
- `td_set` [`ArrheniusTD`](@ref) or [`ArrheniusPeakTD`](@ref) type struct
- `T` Leaf temperature in `[K]`

The equation used for [`ArrheniusTD`](@ref) is
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
```

The equations used for [`ArrheniusPeakTD`](@ref) are
```math
corr = \\exp \\left( \\dfrac{ΔHa}{R T_0} - \\dfrac{ΔHa}{R T_1} \\right)
       \\cdot
       \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
              { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```

The equation used for [`Q10TD`](@ref) is
```math
corr = \\left( \\dfrac{T_1 - T_\\text{REF}}{10} \\right)^{Q_{10}}
```
"""
function temperature_correction(
            td_set::Arrhenius{FT},
            T::FT
) where {FT<:AbstractFloat}
    @unpack ΔHA, T_REF = td_set;

    return exp( ΔHA / GAS_R(FT) * (1 / T_REF - 1 / T) )
end




function temperature_correction(
            td_set::ArrheniusPeak{FT},
            T::FT
) where {FT<:AbstractFloat}
    @unpack ΔHA, ΔHD, ΔSV, T_REF = td_set;

    # _f_a: activation correction, C/_f_b: de-activation correction
    _f_a::FT = exp( ΔHA / GAS_R(FT) * (1 / T_REF - 1 / T) );
    _f_b::FT = (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * T_REF))) / (1 + exp(ΔSV / GAS_R(FT) - ΔHD / (GAS_R(FT) * T)));

    return _f_a * _f_b
end




function temperature_correction(
            td_set::Q10{FT},
            T::FT
) where {FT<:AbstractFloat}
    @unpack Q_10, T_REF = td_set;

    return Q_10 ^ ( (T - T_REF) / 10 )
end
