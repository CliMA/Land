#=
The structure tree of AbstractKpepTD
AbstractKpepTD
---> struct KpepTDBoyd    # using Arrhenius correction
---> struct KpepTDCLM     # using Arrhenius correction
=#
abstract type AbstractKpepTD end




"""
    KpepTDBoyd{FT} <: AbstractKpepTD

A non-mutable AbstractKpepTD type `KpepTDBoyd` that stores information for Kpep temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Boyd et al. (2001) "Temperature responses of C4 photosynthesis: biochemical analysis of Rubisco, phosphoenolpyruvate carboxylase, and carbonic anhydrase in Setaria viridis."

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KpepTDBoyd{FT} <: AbstractKpepTD
    "Kpep at 298.15 K `[Pa]`"
    Kpep       ::FT = FT( 16.0 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 36300.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 36300.0 / (GAS_R*K_25) )
end




"""
    KpepTDCLM{FT} <: AbstractKpepTD

A non-mutable AbstractKpepTD type `KpepTDCLM` that stores information for Kpep temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from ""

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct KpepTDCLM{FT} <: AbstractKpepTD
    "Kpep at 298.15 K `[Pa]`"
    Kpep       ::FT = FT( 8.0 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT( 36000.0 / GAS_R )
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT( 36000.0 / (GAS_R*K_25) )
end
