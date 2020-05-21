#=
The structure tree of AbstractKc
AbstractKc
---> struct KcBernacchi    # using Arrhenius correction
---> struct KcCLM          # using Arrhenius correction
=#
abstract type AbstractKc end




"""
    KcBernacchi{FT} <: AbstractKc

A non-mutable AbstractKc type `KcBernacchi` that stores information for Kc temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KcBernacchi{FT} <: AbstractKc
    "Kc at 298.15 K (404.9 ppm) `[Pa]`"
    Kc         ::FT = FT( 41.0264925 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 79430.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 79430.0 / (GAS_R*K_25) )
end




"""
    KcCLM{FT} <: AbstractKc

A non-mutable AbstractKc type `KcCLM` that stores information for Kc temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KcCLM{FT} <: AbstractKc
    "Kc at 298.15 K `[Pa]`"
    Kc         ::FT = FT( 40.49 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 79430.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 79430.0 / (GAS_R*K_25) )
end
