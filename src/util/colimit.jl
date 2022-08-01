#######################################################################################################################################################################################################
#
# Changes to this type
# General
#     2022-Jan-14: add minimum colimitation struct
#     2022-Jan-14: add quadratic colimitation struct
#     2022-Feb-28: add serial colimitation struct
#     2022-Jan-18: add abstract colimitation type
#     2022-Jul-13: add shortcut constructors for QuadraticColimit
#     2022-Aug-01: add square colimitation struct
#
#######################################################################################################################################################################################################
"""

$(TYPEDEF)

Hierarchy of `AbstractColimit`
- [`MinimumColimit`](@ref)
- [`QuadraticColimit`](@ref)
- [`SerialColimit`](@ref)

"""
abstract type AbstractColimit{FT<:AbstractFloat} end


""" Empty structure to indicate minimum colimitation: `x = min(x₁, x₂)` """
struct MinimumColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end


""" Structure to indicate quadratic colimitation (contains field `CURVATURE`): `θ⋅x² - (x₁ + x₂)⋅x + x₁x₂ = 0` """
Base.@kwdef mutable struct QuadraticColimit{FT<:AbstractFloat} <: AbstractColimit{FT}
    "Curvature factor"
    CURVATURE::FT = 0.98
end


""" Empty structure to indicate serial colimitation: `x = 1 / (1/x₁ + 1/x₂)` """
struct SerialColimit{FT<:AbstractFloat} <:AbstractColimit{FT} end


""" Empty structure to indicate square colimitation: `x = x₁⋅x₂ / sqrt(x₁² + x₂²)` """
struct SquareColimit{FT<:AbstractFloat} <: AbstractColimit{FT} end


ColimitCJCLMC3(FT) = QuadraticColimit{FT}(CURVATURE = 0.98);
ColimitCJCLMC4(FT) = QuadraticColimit{FT}(CURVATURE = 0.8);
ColimitIPCLM(FT) = QuadraticColimit{FT}(CURVATURE = 0.95);
ColimitJCLM(FT)  = QuadraticColimit{FT}(CURVATURE = 0.7);
