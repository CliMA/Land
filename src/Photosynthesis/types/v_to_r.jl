#=
The structure tree of AbstractVtoR
AbstractVtoR
---> VtoRCollatz    # A ratio of 0.01
=#
abstract type AbstractVtoR end




"""
    VtoRCollatz{FT} <: AbstractVtoR

A non-mutable AbstractVtoR type `VtoRCollatz` that stores information for Respiration Vcmax25 correction.
The equation used is `r25 = v25 * correction`.
The data source for the constants can be referred from Collatz et al. (1991) "Physiological and environmental regulation of stomatal conductance, photosynthesis and transpiration: a model that includes a laminar boundary layer."

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VtoRCollatz{FT} <: AbstractVtoR
    "Ratio between r25 and v25"
    v_to_r::FT = FT(0.01)
end