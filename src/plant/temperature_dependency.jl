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
Y_1 = Y_0 \\left( \\dfrac{T_1 - T_0}{10} \\right)^{Q_{10}}
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
