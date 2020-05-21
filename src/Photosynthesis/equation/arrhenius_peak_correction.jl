"""
arrhenius_correction(para_set::JmaxBernacchi{FT}, t_leaf::FT)

Arrhenius correction over Jmax, given
- `para_set` A JmaxBernacchi type parameter set
- `t_leaf` Leaf temperature

The equation used for Arrhenius correction is `f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]`, `f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]`, `C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]`, and `correction = f_above /f_below`
The Arrhenius correction is used for Vcmax and Jmax.
See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model."
"""
function arrhenius_correction(para_set::JmaxBernacchi{FT}, t_leaf::FT) where {FT}
@unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = para_set
return C * exp( ΔHa_to_RT25 * (1-FT(K_25)/t_leaf) ) / ( 1 + exp(ΔSv_to_R - ΔHd_to_R/t_leaf) )
end




"""
    arrhenius_correction(para_set::JmaxCLM{FT}, t_leaf::FT)

Arrhenius correction over Jmax, given
- `para_set` A JmaxCLM type parameter set
- `t_leaf` Leaf temperature

The equation used for Arrhenius correction is `f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]`, `f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]`, `C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]`, and `correction = f_above /f_below`
The Arrhenius correction is used for Vcmax and Jmax.
See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model."
"""
function arrhenius_correction(para_set::JmaxCLM{FT}, t_leaf::FT) where {FT}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = para_set
    return C * exp( ΔHa_to_RT25 * (1-FT(K_25)/t_leaf) ) / ( 1 + exp(ΔSv_to_R - ΔHd_to_R/t_leaf) )
end




"""
    arrhenius_correction(para_set::JmaxLeuning{FT}, t_leaf::FT)

Arrhenius correction over Jmax, given
- `para_set` A JmaxLeuning type parameter set
- `t_leaf` Leaf temperature

The equation used for Arrhenius correction is `f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]`, `f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]`, `C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]`, and `correction = f_above /f_below`
The Arrhenius correction is used for Vcmax and Jmax.
See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model."
"""
function arrhenius_correction(para_set::JmaxLeuning{FT}, t_leaf::FT) where {FT}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = para_set
    return C * exp( ΔHa_to_RT25 * (1-FT(K_25)/t_leaf) ) / ( 1 + exp(ΔSv_to_R - ΔHd_to_R/t_leaf) )
end




"""
    arrhenius_correction(para_set::VcmaxCLM{FT}, t_leaf::FT)

Arrhenius correction over Vcmax, given
- `para_set` A VcmaxCLM type parameter set
- `t_leaf` Leaf temperature

The equation used for Arrhenius correction is `f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]`, `f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]`, `C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]`, and `correction = f_above /f_below`
The Arrhenius correction is used for Vcmax and Jmax.
See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model."
"""
function arrhenius_correction(para_set::VcmaxCLM{FT}, t_leaf::FT) where {FT}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = para_set
    return C * exp( ΔHa_to_RT25 * (1-FT(K_25)/t_leaf) ) / ( 1 + exp(ΔSv_to_R - ΔHd_to_R/t_leaf) )
end




"""
    arrhenius_correction(para_set::VcmaxLeuning{FT}, t_leaf::FT)

Arrhenius correction over Vcmax, given
- `para_set` A VcmaxLeuning type parameter set
- `t_leaf` Leaf temperature

The equation used for Arrhenius correction is `f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]`, `f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]`, `C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]`, and `correction = f_above /f_below`
The Arrhenius correction is used for Vcmax and Jmax.
See Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model."
"""
function arrhenius_correction(para_set::VcmaxLeuning{FT}, t_leaf::FT) where {FT}
    @unpack C, ΔHa_to_RT25, ΔHd_to_R, ΔSv_to_R = para_set
    return C * exp( ΔHa_to_RT25 * (1-FT(K_25)/t_leaf) ) / ( 1 + exp(ΔSv_to_R - ΔHd_to_R/t_leaf) )
end
