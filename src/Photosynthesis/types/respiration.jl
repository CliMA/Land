#=
The structure tree of AbstractRespirationTD
AbstractRespirationTD
---> struct RespirationTDBernacchi    # using Arrhenius correction
---> struct RespirationTDCLM          # using Arrhenius peak correction
=#
abstract type AbstractRespirationTD end




"""
    RespirationTDBernacchi{FT} <: AbstractRespirationTD

A non-mutable AbstractRespirationTD type `RespirationTDBernacchi` that stores information for Respiration temperature correction.
The equation used for temperature correction is `correction = exp( c - ΔHa/(R*T_leaf) )`.
The data source for the constants can be referred from Bernacchi et al. (2001) "Improved temperature response functions for models of Rubisco‐limited photosynthesis".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RespirationTDBernacchi{FT} <: AbstractRespirationTD
    "Ratio between ΔHa and R `[K]`"
    ΔHa_to_R   ::FT = FT(46390.0) / FT(GAS_R)
    "Ratio between ΔHa and R*K-25"
    ΔHa_to_RT25::FT = FT(46390.0) / FT(GAS_R*K_25)
end




"""
    RespirationTDCLM{FT} <: AbstractRespirationTD

A non-mutable AbstractRespirationTD type `RespirationTDCLM` that stores information for Respiration temperature correction.
The equation used for temperature correction is `Arrhenius equation` (see Leuning (2002) "Leuning, R. "Temperature dependence of two parameters in a photosynthesis model").
The data source for the constants can be referred from "".

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RespirationTDCLM{FT} <: AbstractRespirationTD
    "Ratio between ΔHa and R*K_25"
    ΔHa_to_RT25::FT = FT(46390.0) / FT(GAS_R*K_25)
    "Ratio between ΔHd and R"
    ΔHd_to_R   ::FT = FT(150650.0) / FT(GAS_R)
    "Ratio between ΔSv and R"
    ΔSv_to_R   ::FT = FT(490.0) / FT(GAS_R)
    "Correction factor C = 1 + exp( Sv/R + Hd/(RT0) )"
    C          ::FT = 1 + exp( ΔSv_to_R - ΔHd_to_R/FT(K_25) )
end
