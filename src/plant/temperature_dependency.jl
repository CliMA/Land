#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-13: migrate abstract temperature dependency type from Photosynthesis.jl
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of AbstractTemperatureDependency:
- [`Arrhenius`](@ref)
- [`ArrheniusPeak`](@ref)
- [`Q10`](@ref)

"""
abstract type AbstractTemperatureDependency{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to Arrhenius
#     2022-Jan-13: define the struct mutable, use ΔHA directly in the struct, add field T_REF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An [`Arrhenius`](@ref) type struct using
```math
Y_1 = Y_0 \\cdot \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right)
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Arrhenius{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to ArrheniusPeak
#     2022-Jan-13: define the struct mutable, use ΔHA/ΔHD/ΔSV directly in the struct, add field T_REF/VAL_REF
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An [`ArrheniusPeak`](@ref) type struct using
```math
Y_1 = Y_0 \\cdot \\exp \\left( \\dfrac{H_a}{R T_0} - \\dfrac{H_a}{R T_1} \\right)
          \\cdot \\dfrac{ 1 + \\exp \\left( \\dfrac{S_v T_0 - H_d}{R T_0} \\right) }
                        { 1 + \\exp \\left( \\dfrac{S_v T_1 - H_d}{R T_1} \\right) }
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct ArrheniusPeak{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
    "Deactivation energy"
    ΔHD::FT
    "Entropy factor"
    ΔSV::FT
end


#######################################################################################################################################################################################################
#
# Changes to this struct
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to Q10
#     2022-Jan-14: make structure mutable
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

An [`Q10`](@ref) type struct using
```math
Y_1 = Y_0 \\cdot Q_{10} ^ \\dfrac{T_1 - T_0}{10}
```

# Fields

$(TYPEDFIELDS)

"""
Base.@kwdef mutable struct Q10{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    # General model information
    "Power of Q10 correction"
    Q_10::FT
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
end


#######################################################################################################################################################################################################
#
# Changes to the constructors
# General
#     2022-Jan-14: migrate from Photosynthesis.jl
#     2022-Feb-11: add temperature dependent Jmax and Vcmax TD from CLM
# Sources
#     Lavigne and Ryan (1997) Growth and maintenance respiration rates of aspen, blackspruce and jack pine stems at northern and southern BOREAS sites
#     Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco‐limited photosynthesis
#     Boyd et al. (2001) Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
#     Leuning (2002) Temperature dependence of two parameters in a photosynthesis model
#     Kattge et al. (2007) Temperature acclimation in a biochemical model of photosynthesis: a reanalysis of data from 36 species
#     Sperry et al. (2019) The impact of rising CO2 and acclimation on the response of US forests to global warming
#     Johnson et al. (2021) The limiting factors and regulatory processes that control the environmental responses of C3, C3–C4 intermediate, and C4 photosynthesis
#     CLM5 Documentation. Chapter 9 Page 106
#
#######################################################################################################################################################################################################
KcTDBernacchi(FT)          = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 41.0264925, ΔHA = 79430.0);
KcTDCLM(FT)                = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 40.49     , ΔHA = 79430.0);
KoTDBernacchi(FT)          = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 28208.88  , ΔHA = 36380.0);
KoTDCLM(FT)                = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 27840.0   , ΔHA = 36380.0);
KpepTDCLM(FT)              = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 8.0       , ΔHA = 36000.0);
KpepTDBoyd(FT)             = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 16.0      , ΔHA = 36300.0);
KqTDJohnson(FT)            = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 300       , ΔHA = 37000.0);
RespirationTDBernacchi(FT) = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = NaN       , ΔHA = 46390.0);
VcmaxTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = NaN       , ΔHA = 65330.0);
VomaxTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = NaN       , ΔHA = 60110.0);
ΓStarTDBernacchi(FT)       = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 4.33164375, ΔHA = 37830.0);
ΓStarTDCLM(FT)             = Arrhenius{FT}(T_REF = T₂₅(), VAL_REF = 4.275     , ΔHA = 37830.0);

JmaxTDBernacchi(FT)               = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 57500.0, ΔHD = 439000.0, ΔSV = 1400.0);
JmaxTDCLM(FT, t::Number = T₂₅())  = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 50000.0, ΔHD = 200000.0, ΔSV = 659.70 - 0.75 * (t - T₀()) );
JmaxTDLeuning(FT)                 = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 50300.0, ΔHD = 152044.0, ΔSV = 495.0 );
RespirationTDCLM(FT)              = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 46390.0, ΔHD = 150650.0, ΔSV = 490.0 );
VcmaxTDCLM(FT, t::Number = T₂₅()) = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 72000.0, ΔHD = 200000.0, ΔSV = 668.39 - 1.07 * (t - T₀()) );
VcmaxTDLeuning(FT)                = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 73637.0, ΔHD = 149252.0, ΔSV = 486.0 );
VpmaxTDBoyd(FT)                   = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = NaN , ΔHA = 94800.0, ΔHD = 73300.0 , ΔSV = 250.0 );
ΗCTDJohnson(FT)                   = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = 1.0 , ΔHA = 0.0    , ΔHD = 220000.0, ΔSV = 710.0 );
ΗLTDJohnson(FT)                   = ArrheniusPeak{FT}(T_REF = T₂₅(), VAL_REF = 0.75, ΔHA = 0.0    , ΔHD = 220000.0, ΔSV = 710.0 );

Q10TDAngiosperm(FT) = Q10{FT}(Q_10 = 1.4, T_REF = T₂₅(), VAL_REF = 0.0140/8760);
Q10TDGymnosperm(FT) = Q10{FT}(Q_10 = 1.7, T_REF = T₂₅(), VAL_REF = 0.0425/8760);
