#=
The structure tree of AbstractΓStar
AbstractΓStar
---> struct ΓStarBernacchi    # using Arrhenius correction
---> struct ΓStarCLM          # using Arrhenius correction
=#
abstract type AbstractΓStar end




"""
    ΓStarBernacchi{FT} <: AbstractΓStar

A non-mutable AbstractΓStar type `ΓStarBernacchi` that stores information for Γ_star temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ΓStarBernacchi{FT} <: AbstractΓStar
    "Γ_star at 298.15 K (42.75 ppm) `[Pa]`"
    Γ_star     ::FT = FT(4.33164375)
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(37830.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(37830.0) / FT(GAS_R*K_25)
end




"""
    ΓStarCLM{FT} <: AbstractΓStar

A non-mutable AbstractΓStar type `ΓStarCLM` that stores information for Γ_star temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ΓStarCLM{FT} <: AbstractΓStar
    "Γ_star at 298.15 K `[Pa]`"
    Γ_star     ::FT = FT( 4.275 )
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(37830.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(37830.0) / FT(GAS_R*K_25)
end
