#=
The structure tree of AbstractKcTD
AbstractKcTD
---> struct KcTDBernacchi    # using Arrhenius correction
---> struct KcTDCLM          # using Arrhenius correction
=#
abstract type AbstractKcTD end




"""
    KcTDBernacchi{FT} <: AbstractKcTD

A non-mutable AbstractKcTD type `KcTDBernacchi` that stores information for Kc temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KcTDBernacchi{FT} <: AbstractKcTD
    "Kc at 298.15 K (404.9 ppm) `[Pa]`"
    Kc         ::FT = FT( 41.0264925 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 79430.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 79430.0 / (GAS_R*K_25) )
end




"""
    KcTDCLM{FT} <: AbstractKcTD

A non-mutable AbstractKcTD type `KcTDCLM` that stores information for Kc temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KcTDCLM{FT} <: AbstractKcTD
    "Kc at 298.15 K `[Pa]`"
    Kc         ::FT = FT( 40.49 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 79430.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 79430.0 / (GAS_R*K_25) )
end
