"""
$(TYPEDEF)

Hierachy of AbstractTemperatureDependency:
- [`Arrhenius`](@ref)
- [`ArrheniusPeak`](@ref)
"""
abstract type AbstractTemperatureDependency{FT<:AbstractFloat} end


#######################################################################################################################################################################################################
#
# Changes to this structure
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
mutable struct Arrhenius{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
end


#######################################################################################################################################################################################################
#
# Changes to this structure
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
mutable struct ArrheniusPeak{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Activation energy"
    ΔHA::FT
    "Deactivation energy"
    ΔHD::FT
    "Entrophy factor"
    ΔSV::FT
end


#######################################################################################################################################################################################################
#
# Changes to this structure
# General
#     2022-Jan-13: migrate from Photosynthesis.jl, rename to Q10
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
struct Q10{FT<:AbstractFloat} <: AbstractTemperatureDependency{FT}
    "Reference temperature `[K]`"
    T_REF::FT
    "Uncorrected vakye at reference temperature"
    VAL_REF::FT
    "Power of Q10 correction"
    Q_10::FT
end


#######################################################################################################################################################################################################
#
# Changes to the constructors
# General
#     2022-Jan-14: migrate from Photosynthesis.jl
# Sources
#     Lavigne and Ryan (1997) Growth and maintenance respiration rates of aspen, blackspruce and jack pine stems at northern and southern BOREAS sites
#     Bernacchi et al. (2001) Improved temperature response functions for models of Rubisco‐limited photosynthesis
#     Boyd et al. (2001) Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis
#     Leuning (2002) Temperature dependence of two parameters in a photosynthesis model
#
#######################################################################################################################################################################################################
KcTDBernacchi(FT)          = Arrhenius{FT}(T_25(), 41.0264925, 79430.0);
KcTDCLM(FT)                = Arrhenius{FT}(T_25(), 40.49     , 79430.0);
KoTDBernacchi(FT)          = Arrhenius{FT}(T_25(), 28208.88  , 36380.0);
KoTDCLM(FT)                = Arrhenius{FT}(T_25(), 27840.0   , 36380.0);
KpepTDCLM(FT)              = Arrhenius{FT}(T_25(), 8.0       , 36000.0);
KpepTDBoyd(FT)             = Arrhenius{FT}(T_25(), 16.0      , 36300.0);
RespirationTDBernacchi(FT) = Arrhenius{FT}(T_25(), NaN       , 46390.0);
VcmaxTDBernacchi(FT)       = Arrhenius{FT}(T_25(), NaN       , 65330.0);
VomaxTDBernacchi(FT)       = Arrhenius{FT}(T_25(), NaN       , 60110.0);
ΓStarTDBernacchi(FT)       = Arrhenius{FT}(T_25(), 4.33164375, 37830.0);
ΓStarTDCLM(FT)             = Arrhenius{FT}(T_25(), 4.275     , 37830.0);

JmaxTDBernacchi(FT)  = ArrheniusPeak{FT}(T_25(), NaN, 57500.0, 439000.0, 1400.0);
JmaxTDCLM(FT)        = ArrheniusPeak{FT}(T_25(), NaN, 43540.0, 150000.0, 490.0 );
JmaxTDLeuning(FT)    = ArrheniusPeak{FT}(T_25(), NaN, 50300.0, 152044.0, 495.0 );
RespirationTDCLM(FT) = ArrheniusPeak{FT}(T_25(), NaN, 46390.0, 150650.0, 490.0 );
VcmaxTDCLM(FT)       = ArrheniusPeak{FT}(T_25(), NaN, 65330.0, 150000.0, 490.0 );
VcmaxTDLeuning(FT)   = ArrheniusPeak{FT}(T_25(), NaN, 73637.0, 149252.0, 486.0 );
VpmaxTDBoyd(FT)      = ArrheniusPeak{FT}(T_25(), NaN, 94800.0, 73300.0 , 250.0 );

Q10TDAngiosperm(FT) = Q10{FT}(T_25(), 0.0140/8760, 1.4);
Q10TDGymnosperm(FT) = Q10{FT}(T_25(), 0.0425/8760, 1.7);
