"""
    arrhenius_correction(para_set::KcBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KcBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."
"""
function arrhenius_correction(para_set::KcBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KcCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KcCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
"""
function arrhenius_correction(para_set::KcCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KoBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KoBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."
"""
function arrhenius_correction(para_set::KoBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KoCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KoCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
"""
function arrhenius_correction(para_set::KoCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::RespirationBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `RespirationBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."
"""
function arrhenius_correction(para_set::RespirationBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::VcmaxBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `VcmaxBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."
"""
function arrhenius_correction(para_set::VcmaxBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::VomaxBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `VcomaxBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."
"""
function arrhenius_correction(para_set::VomaxBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::ΓStarBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `ΓStarBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."
"""
function arrhenius_correction(para_set::ΓStarBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::ΓStarCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `ΓStarCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Rd, and Γ_star.
"""
function arrhenius_correction(para_set::ΓStarCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end
