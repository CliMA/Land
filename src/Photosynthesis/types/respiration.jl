#=
The structure tree of AbstractRespiration
AbstractRespiration
---> struct RespirationBernacchi    # using Arrhenius correction
=#
abstract type AbstractRespiration end




"""
    RespirationBernacchi{FT} <: AbstractRespiration

A non-mutable AbstractRespiration type `RespirationBernacchi` that stores information for Respiration temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RespirationBernacchi{FT} <: AbstractRespiration
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(46390.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(46390.0) / FT(GAS_R*K_25)
end
