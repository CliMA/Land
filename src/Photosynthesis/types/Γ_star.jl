#=
The structure tree of AbstractΓStarTD
AbstractΓStarTD
---> struct ΓStarTDBernacchi    # using Arrhenius correction
---> struct ΓStarTDCLM          # using Arrhenius correction
=#
abstract type AbstractΓStarTD end




"""
    ΓStarTDBernacchi{FT} <: AbstractΓStarTD

A non-mutable AbstractΓStarTD type `ΓStarTDBernacchi` that stores information for Γ_star temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ΓStarTDBernacchi{FT} <: AbstractΓStarTD
    "Γ_star at 298.15 K (42.75 ppm) `[Pa]`"
    Γ_star     ::FT = FT(4.33164375)
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(37830.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(37830.0) / FT(GAS_R*K_25)
end




"""
    ΓStarTDCLM{FT} <: AbstractΓStarTD

A non-mutable AbstractΓStarTD type `ΓStarTDCLM` that stores information for Γ_star temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ΓStarTDCLM{FT} <: AbstractΓStarTD
    "Γ_star at 298.15 K `[Pa]`"
    Γ_star     ::FT = FT( 4.275 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(37830.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(37830.0) / FT(GAS_R*K_25)
end
