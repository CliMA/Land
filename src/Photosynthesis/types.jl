###############################################################################
#
# Temperature dependency parameter set
#
###############################################################################
"""
    AbstractTDParameterSet

Hierarchy of the `AbstractTDParameterSet`:
- AbstractArrheniusTDParameterSet
- AbstractArrheniusPeakTDParameterSet
"""
abstract type AbstractTDParameterSet end




"""
    AbstractArrheniusTDParameterSet

A parameter set type that uses Arrhenius correction using `correction = exp( c - ΔHa/(RT) )`.
"""
abstract type AbstractArrheniusTDParameterSet <: AbstractTDParameterSet end




"""
    struct ArrheniusTD{FT,val_25,ΔHa}
    
An `AbstractArrheniusTDParameterSet` type struct that takes the following to initialize
- `FT` Floating type
- `val_25` Uncorrected value at 298.15 K
- `ΔHa` Heat of corrosion reaction

# Fields
$(DocStringExtensions.FIELDS)
"""
@Base.kwdef struct ArrheniusTD{FT,val_25,ΔHa} <: AbstractArrheniusTDParameterSet
    "Uncorrected value at 298.15 K"
    VAL_25::FT = FT( val_25 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R::FT = FT( ΔHa / GAS_R )
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT( ΔHa / (GAS_R*K_25) )
end




"""
    AbstractArrheniusPeakTDParameterSet

A parameter set type that uses Arrhenius peak correction using
```
    C = 1 + exp[(Sv*T_0 - Hd )/(RT_0 )]
    f_above = C * exp [ ( Ha /RT_0 ) * ( 1 - T_0/T_1 ) ]
    f_below = 1 + exp [ ( Sv*T_1 - Hd ) / ( RT_1) ]
    correction = f_above /f_below
```
"""
abstract type AbstractArrheniusPeakTDParameterSet <: AbstractTDParameterSet end




"""
    struct ArrheniusPeakTD{FT,val_25,ΔHa}
    
An `AbstractArrheniusPeakTDParameterSet` type struct that takes the following to initialize
- `FT` Floating type
- `ΔHa` Heat of corrosion reaction
- `ΔHd` Dimerization enthalpy
- `ΔSv` Entropy

# Fields
$(DocStringExtensions.FIELDS)
"""
@Base.kwdef struct ArrheniusPeakTD{FT,ΔHa,ΔHd,ΔSv} <: AbstractArrheniusPeakTDParameterSet
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT( ΔHa / (GAS_R*K_25) )
    "Ratio between ΔHd and R"
    ΔHd_to_R::FT = FT( ΔHd / GAS_R )
    "Ratio between ΔSv and R"
    ΔSv_to_R::FT = FT( ΔSv / GAS_R )
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end








###############################################################################
#
# Photosynthesis model parameter set -- temperature dependencies and etc
#
###############################################################################
"""
    AbstractPhotoModelParaSet

A parameter set that stores TD parasets and constants
"""
abstract type AbstractPhotoModelParaSet end




"""
    struct C3Paraset{FT, JTD, KcTD, KoTD, ReTD, VcTD, ΓsTD, V2R, eff_1, eff_2}

An `AbstractPhotoModelParaSet` type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct C3ParaSet{FT, JTD, KcTD, KoTD, ReTD, VcTD, ΓsTD, V2R, eff_1, eff_2} <: AbstractPhotoModelParaSet
    "Jmax temperature dependency"
    JT ::AbstractTDParameterSet = JTD
    "Kc temperature dependency"
    KcT::AbstractTDParameterSet = KcTD
    "Ko temperature dependency"
    KoT::AbstractTDParameterSet = KoTD
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet = ReTD
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet = VcTD
    "Γ_star temperature dependency"
    ΓsT::AbstractTDParameterSet = ΓsTD
    "Vcmax25 and respiration correlation"
    VR::FT = V2R
    "Coefficient 4.0/4.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_1::FT = eff_1
    "Coefficient 8.0/10.5 for NADPH/ATP requirement stochiometry, respectively"
    Eff_2::FT = eff_2
end




"""
    struct C4ParaSet{FT, KpTD, ReTD, VcTD, VpTD, V2R}

An `AbstractPhotoModelParaSet` type of parameter sets for different temperature dependencies.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct C4ParaSet{FT, KpTD, ReTD, VcTD, VpTD, V2R} <: AbstractPhotoModelParaSet
    "Kpep temperature dependency"
    KpT::AbstractTDParameterSet = KpTD
    "Respiration temperature dependency"
    ReT::AbstractTDParameterSet = ReTD
    "Vcmax temperature dependency"
    VcT::AbstractTDParameterSet = VcTD
    "Vpmax temperature dependency"
    VpT::AbstractTDParameterSet = VpTD
    "Vcmax25 and respiration correlation"
    VR::FT = V2R
end








###############################################################################
#
# Leaf fluorescence-related parameter set
#
###############################################################################
"""
    AbstractFluoModelParaSet

A parameter set that stores fluorescence model parameters
"""
abstract type AbstractFluoModelParaSet end




"""
    struct FluoParaSet{FT, kn1, kn2, kn3}

A `AbstractFluoModelParaSet` type paramter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct FluoParaSet{FT, kn1, kn2, kn3} <: AbstractFluoModelParaSet
    "Fluorescence model coefficient"
    Kn1::FT = kn1
    "Fluorescence model coefficient"
    Kn2::FT = kn2
    "Fluorescence model coefficient"
    Kn3::FT = kn3
end
