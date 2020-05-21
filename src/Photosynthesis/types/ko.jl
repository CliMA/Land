#=
The structure tree of AbstractKoTD
AbstractKoTD
---> struct KoTDBernacchi    # using Arrhenius correction
---> struct KoTDCLM          # using Arrhenius correction
=#
abstract type AbstractKoTD end




"""
    KoBernacchi{FT} <: AbstractKoTD

A non-mutable AbstractKoTD type `KoTDBernacchi` that stores information for Ko temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KoTDBernacchi{FT} <: AbstractKoTD
    "Ko at 298.15 K (278.4 ‰) `[Pa]`"
    Ko         ::FT = FT( 28208.88 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 36380.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 36380.0 / (GAS_R*K_25) )
end




"""
    KoTDCLM{FT} <: AbstractKoTD

A non-mutable AbstractKoTD type `KoTDCLM` that stores information for Ko temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KoTDCLM{FT} <: AbstractKoTD
    "Ko at 298.15 K `[Pa]`"
    Ko         ::FT = FT( 27840.0 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 36380.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 36380.0 / (GAS_R*K_25) )
end

