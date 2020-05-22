"""
    arrhenius_correction(para_set::KcTDBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KcTDBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."

"""
function arrhenius_correction(para_set::KcTDBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KcTDCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KcTDCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.

"""
function arrhenius_correction(para_set::KcTDCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KoTDBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KoTDBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."

"""
function arrhenius_correction(para_set::KoTDBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KoTDCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KoTDCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.

"""
function arrhenius_correction(para_set::KoTDCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KpepTDBoyd{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KpepTDBoyd` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.

"""
function arrhenius_correction(para_set::KpepTDBoyd{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::KpepTDCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `KpepTDCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.

"""
function arrhenius_correction(para_set::KpepTDCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::RespirationTDBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `RespirationTDBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."

"""
function arrhenius_correction(para_set::RespirationTDBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::VcmaxTDBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `VcmaxTDBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, kpep, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."

"""
function arrhenius_correction(para_set::VcmaxTDBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::VomaxTDBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `VomaxTDBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."

"""
function arrhenius_correction(para_set::VomaxTDBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::ΓStarTDBernacchi{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `ΓStarTDBernacchi` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.
See Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco-limited photosynthesis."

"""
function arrhenius_correction(para_set::ΓStarTDBernacchi{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end




"""
    arrhenius_correction(para_set::ΓStarTDCLM{FT}, t_leaf::FT)

A correction factor based on arrhenius's fitting procedure, given
- `para_set` A `ΓStarTDCLM` type parameter set
- `t_leaf` Leaf temperature in `[K]`

The equation used is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The arrhenius correction can be used for Vcmax, Vomax, Kc, Ko, Kpep, Rd, and Γ_star.
    
"""
function arrhenius_correction(para_set::ΓStarTDCLM{FT}, t_leaf::FT) where {FT}
    return exp( -para_set.ΔHa_to_R/t_leaf + para_set.ΔHa_to_RT25 )
end
