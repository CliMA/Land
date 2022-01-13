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
            td_set::ArrheniusTD{FT},
            T::FT
) where {FT<:AbstractFloat}
    return exp( td_set.ΔHa_to_RT25 - td_set.ΔHa_to_R/T )
end




function temperature_correction(
            td_set::ArrheniusPeakTD{FT},
            T::FT
) where {FT<:AbstractFloat}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = td_set;

    # _f_a: activation correction, C/_f_b: de-activation correction
    _f_a::FT = exp( ΔHa_to_RT25 * (1 - T_25(FT)/T) );
    _f_b::FT = 1 + exp(ΔSv_to_R - ΔHd_to_R/T);

    return C /_f_b * _f_a
end




function temperature_correction(
            td_set::Q10TD{FT},
            T::FT
) where {FT<:AbstractFloat}
    return (td_set.Q_10) ^ ( (T - td_set.T_REF) / 10 )
end
